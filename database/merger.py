import sqlite3
import pandas as pd
import os
from pathlib import Path
from config import OUTPUTS_DIR, PROJECT_ROOT


class DatabaseMerger:
    def __init__(self, output_db_path, sql_seed_path=None):
        self.output_db_path = output_db_path

        if sql_seed_path:
            self.sql_seed = open(sql_seed_path, 'r').read()
            self.create_output_db()

        self.conn = sqlite3.connect(self.output_db_path)
        self.cursor = self.conn.cursor()

        self.node_id_mapping = {}

        self.processed_nodes = set()

        self.cursor.execute("PRAGMA table_info(node)")
        self.node_columns = [col[1] for col in self.cursor.fetchall()]
        print(f"Node table columns: {self.node_columns}")

        self.cursor.execute("PRAGMA table_info(edge)")
        self.edge_columns = [col[1] for col in self.cursor.fetchall()]
        print(f"Edge table columns: {self.edge_columns}")

    def create_output_db(self):
        conn = sqlite3.connect(self.output_db_path)
        conn.executescript(self.sql_seed)
        conn.commit()
        conn.close()

    def merge_database(self, db_path):
        print(f"Processing database: {db_path}")

        source_conn = sqlite3.connect(db_path)
        source_conn.row_factory = sqlite3.Row

        self.import_nodes(source_conn)

        self.import_node_identifiers(source_conn)

        self.import_edges(source_conn)

        source_conn.close()
        print(f"Finished processing database: {db_path}")

    def get_node_signature(self, node_data, node_identifiers=None):
        signature_parts = []

        if node_data['name']:
            signature_parts.append(f"name:{node_data['name']}")

        if node_identifiers:
            for id_type, id_value in node_identifiers.items():
                if id_value:
                    signature_parts.append(f"{id_type}:{id_value}")

        return "|".join(sorted(signature_parts))

    def get_node_identifiers_from_db(self, source_conn, node_id):
        cursor = source_conn.cursor()
        cursor.execute("SELECT id_type, id_value FROM node_identifier WHERE node_id = ?", (node_id,))
        return {row['id_type']: row['id_value'] for row in cursor.fetchall()}

    def find_existing_node(self, node_data, node_identifiers):
        self.cursor.execute("SELECT id FROM node WHERE name = ?", (node_data['name'],))
        result = self.cursor.fetchone()
        if result:
            return result[0]

        for id_type, id_value in node_identifiers.items():
            self.cursor.execute("""
                SELECT n.id
                FROM node n
                JOIN node_identifier ni ON n.id = ni.node_id
                WHERE ni.id_type = ? AND ni.id_value = ?
            """, (id_type, id_value))
            result = self.cursor.fetchone()
            if result:
                return result[0]

        return None

    def node_exists(self, signature):
        return signature in self.processed_nodes

    def merge_node_properties(self, existing_node, new_node):
        def merge_strings(str1, str2, separator="|"):
            if not str1 and not str2:
                return ""
            if not str1:
                return str2
            if not str2:
                return str1

            list1 = [item for item in str1.split(separator) if item and item != '-']
            list2 = [item for item in str2.split(separator) if item and item != '-']

            merged = list(set(list1 + list2))
            return separator.join(merged)

        self.cursor.execute("SELECT * FROM node WHERE id = ?", (existing_node,))
        existing_data = self.cursor.fetchone()

        pathways_idx = 6
        source_idx = 7
        function_idx = 8

        merged_pathways = merge_strings(existing_data[pathways_idx], new_node['pathways'])
        merged_source = merge_strings(existing_data[source_idx], new_node['source'])
        merged_function = merge_strings(existing_data[function_idx], new_node['function'])

        self.cursor.execute("""
            UPDATE node
            SET pathways = ?, function = ?, source = ?
            WHERE id = ?
        """, (merged_pathways, merged_function, merged_source, existing_node))

        return existing_node

    def import_nodes(self, source_conn):
        cursor = source_conn.cursor()
        cursor.execute("SELECT * FROM node")
        nodes = cursor.fetchall()

        nodes_added = 0
        nodes_merged = 0

        for node in nodes:
            node_data = {key: node[key] for key in node.keys()}

            identifiers = self.get_node_identifiers_from_db(source_conn, node_data['id'])

            existing_node_id = self.find_existing_node(node_data, identifiers)

            if existing_node_id:
                self.merge_node_properties(existing_node_id, node_data)

                self.node_id_mapping[f"{db_path}_{node_data['id']}"] = existing_node_id

                nodes_merged += 1
            else:
                self.cursor.execute(
                    """
                    INSERT INTO node
                    (name, primary_id_type, display_name, tax_id, type, pathways, source, function)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        node_data['name'],
                        node_data['primary_id_type'],
                        node_data['display_name'],
                        node_data['tax_id'],
                        node_data['type'],
                        node_data['pathways'],
                        node_data['source'],
                        node_data['function']
                    )
                )

                new_node_id = self.cursor.lastrowid

                self.node_id_mapping[f"{db_path}_{node_data['id']}"] = new_node_id

                nodes_added += 1

            if (nodes_added + nodes_merged) % 1000 == 0:
                self.conn.commit()
                print(f"Processed {nodes_added + nodes_merged} nodes so far...")

        self.conn.commit()
        print(f"Imported {nodes_added} new nodes, merged {nodes_merged} existing nodes")

    def import_node_identifiers(self, source_conn):
        """Import node identifiers, mapping to new node IDs."""
        print("Importing node identifiers...")

        cursor = source_conn.cursor()
        cursor.execute("SELECT * FROM node_identifier")
        identifiers = cursor.fetchall()
        print(f"Found {len(identifiers)} node identifiers in source database")

        identifiers_added = 0
        identifiers_skipped = 0

        for identifier in identifiers:
            old_node_id = identifier['node_id']
            id_type = identifier['id_type']
            id_value = identifier['id_value']
            is_primary = identifier['is_primary']

            new_node_id = self.node_id_mapping.get(f"{db_path}_{old_node_id}")

            if new_node_id:
                self.cursor.execute(
                    "SELECT 1 FROM node_identifier WHERE node_id = ? AND id_type = ?",
                    (new_node_id, id_type)
                )
                exists = self.cursor.fetchone()

                if exists:
                    self.cursor.execute(
                        """
                        UPDATE node_identifier
                        SET id_value = ?, is_primary = ?
                        WHERE node_id = ? AND id_type = ?
                        """,
                        (id_value, is_primary, new_node_id, id_type)
                    )
                    identifiers_skipped += 1
                else:
                    try:
                        self.cursor.execute(
                            """
                            INSERT INTO node_identifier
                            (node_id, id_type, id_value, is_primary)
                            VALUES (?, ?, ?, ?)
                            """,
                            (new_node_id, id_type, id_value, is_primary)
                        )
                        identifiers_added += 1
                    except sqlite3.IntegrityError:
                        identifiers_skipped += 1

            if (identifiers_added + identifiers_skipped) % 1000 == 0:
                self.conn.commit()
                print(f"Processed {identifiers_added + identifiers_skipped} identifiers so far...")

        self.conn.commit()
        print(f"Imported {identifiers_added} new identifiers, skipped {identifiers_skipped} existing ones")

    def merge_edge_properties(self, existing_edge_id, new_edge_data):
        def merge_strings(str1, str2, separator="|"):
            if not str1 and not str2:
                return ""
            if not str1:
                return str2
            if not str2:
                return str1

            list1 = [item for item in str1.split(separator) if item and item != '-']
            list2 = [item for item in str2.split(separator) if item and item != '-']

            merged = list(set(list1 + list2))
            return separator.join(merged)

        self.cursor.execute("SELECT * FROM edge WHERE id = ?", (existing_edge_id,))
        existing_data = self.cursor.fetchone()

        interaction_types_idx = 6
        source_db_idx = 7

        merged_interaction_types = merge_strings(existing_data[interaction_types_idx], new_edge_data['interaction_types'])
        merged_source_db = merge_strings(existing_data[source_db_idx], new_edge_data['source_db'])

        self.cursor.execute("""
            UPDATE edge
            SET interaction_types = ?, source_db = ?
            WHERE id = ?
        """, (merged_interaction_types, merged_source_db, existing_edge_id))

    def edge_exists(self, source_id, target_id, layer):
        self.cursor.execute("""
            SELECT id FROM edge
            WHERE interactor_a_node_id = ? AND interactor_b_node_id = ? AND layer = ?
        """, (source_id, target_id, layer))
        result = self.cursor.fetchone()
        if result:
            return result[0]

        self.cursor.execute("""
            SELECT id, interaction_types FROM edge
            WHERE interactor_a_node_id = ? AND interactor_b_node_id = ? AND layer = ?
        """, (target_id, source_id, layer))
        result = self.cursor.fetchone()
        if result:
            interaction_types = result[1]
            if interaction_types and 'is_directed:true' not in interaction_types.split('|'):
                return result[0]

        return None

    def import_edges(self, source_conn):
        print("Importing edges...")

        cursor = source_conn.cursor()
        cursor.execute("SELECT * FROM edge")
        edges = cursor.fetchall()
        print(f"Found {len(edges)} edges in source database")

        edges_added = 0
        edges_merged = 0
        edges_skipped = 0

        for edge in edges:
            old_source_id = edge['interactor_a_node_id']
            old_target_id = edge['interactor_b_node_id']

            new_source_id = self.node_id_mapping.get(f"{db_path}_{old_source_id}")
            new_target_id = self.node_id_mapping.get(f"{db_path}_{old_target_id}")

            if not new_source_id or not new_target_id:
                edges_skipped += 1
                continue

            edge_data = {
                'interactor_a_node_id': new_source_id,
                'interactor_b_node_id': new_target_id,
                'interactor_a_node_name': edge['interactor_a_node_name'],
                'interactor_b_node_name': edge['interactor_b_node_name'],
                'layer': edge['layer'],
                'interaction_types': edge['interaction_types'],
                'source_db': edge['source_db']
            }

            existing_edge_id = self.edge_exists(new_source_id, new_target_id, edge['layer'])

            if existing_edge_id:
                self.merge_edge_properties(existing_edge_id, edge_data)
                edges_merged += 1
            else:
                self.cursor.execute(
                    """
                    INSERT INTO edge
                    (interactor_a_node_id, interactor_b_node_id, interactor_a_node_name,
                    interactor_b_node_name, layer, interaction_types, source_db)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        new_source_id,
                        new_target_id,
                        edge['interactor_a_node_name'],
                        edge['interactor_b_node_name'],
                        edge['layer'],
                        edge['interaction_types'],
                        edge['source_db']
                    )
                )

                edges_added += 1

            if (edges_added + edges_merged + edges_skipped) % 1000 == 0:
                self.conn.commit()
                print(f"Processed {edges_added + edges_merged + edges_skipped} edges so far...")

        self.conn.commit()
        print(f"Imported {edges_added} new edges, merged {edges_merged} existing edges, skipped {edges_skipped} due to missing nodes")

    def finalize(self):
        self.conn.commit()
        self.conn.close()
        print(f"Merged database saved to {self.output_db_path}")


if __name__ == "__main__":
    SQL_SEED = PROJECT_ROOT / "database" / "network_db_seed2.sql"

    OUTPUT_DB = OUTPUTS_DIR / "merged_sources_ok.db"

    input_dbs = [
        OUTPUTS_DIR / "kegg.db",
        OUTPUTS_DIR / "ferrdb_network.db",
        OUTPUTS_DIR / "ferreg_network.db"
    ]

    merger = DatabaseMerger(OUTPUT_DB, SQL_SEED)

    for db_path in input_dbs:
        if os.path.exists(db_path):
            merger.merge_database(db_path)
        else:
            print(f"Warning: Database {db_path} does not exist, skipping")

    merger.finalize()

    print("\nMerge Complete!")
    print(f"Output database: {OUTPUT_DB}")
