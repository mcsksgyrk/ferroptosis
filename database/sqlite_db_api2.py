import sqlite3


class PsimiSQL:
    def __init__(self, sql_seed_file_location):
        self.sql_seed = open(sql_seed_file_location).read()

        self.db = self.create_db(":memory:")
        self.cursor = self.db.cursor()

    def import_from_db_file(self, db_file_location):
        temporary_db = sqlite3.connect(db_file_location)
        temporary_db_name = "file_db"

        params = (db_file_location, temporary_db_name)
        self.db.execute("ATTACH ? AS ? ", params)
        self.db.commit()

        self.db.execute("INSERT INTO node SELECT * FROM %s.node" % temporary_db_name)
        self.db.commit()

        self.db.execute("INSERT INTO edge SELECT * FROM %s.edge" % temporary_db_name)
        self.db.commit()

        try:
            self.db.execute("INSERT INTO node_identifier SELECT * FROM %s.node_identifier" % temporary_db_name)
            self.db.commit()
        except sqlite3.OperationalError:
            pass

        temporary_db.close()

        self.db.execute("DETACH DATABASE %s" % temporary_db_name)
        self.db.commit()

    def create_db(self, location):
        db = sqlite3.connect(location)
        db.text_factory = str

        create_tables_query = self.sql_seed

        db.executescript(create_tables_query)
        db.commit()

        return db

    def check_if_node_exists(self, node_dict):
        query = """
            SELECT n.id FROM node n
            JOIN node_identifier ni ON n.id = ni.node_id
            WHERE ni.id_value = ?
            AND (? IS NULL OR n.tax_id = ?)
            LIMIT 1
            """
        res = self.cursor.execute(query, (node_dict['name'],
                                          node_dict['tax_id'],
                                          node_dict['tax_id'])
                                  ).fetchone()
        return res[0] if res else None

    def check_node_dict_identifiers(self, node_dict):
        for id_type in ['kegg_id', 'uniprot_id', 'pubchem_id', "pubmed_id", "hgnc_id", "ensg_id"]:
            if id_type in node_dict and node_dict[id_type]:
                is_primary = False
                if node_dict['primary_id_type'] == id_type:
                    is_primary = True
                self.insert_node_identifier(
                        node_dict['id'],
                        id_type,
                        node_dict[id_type],
                        is_primary
                )

    def insert_node(self, node_dict):
        node_id = self.check_if_node_exists(node_dict)
        existing_node = None

        if node_id:
            existing_node = self.get_node_by_id(node_id)

        if ('id' not in node_dict) and not existing_node:
            if 'display_name' not in node_dict:
                node_dict['display_name'] = node_dict.get('name', '')

            if 'type' not in node_dict:
                node_dict['type'] = 'protein'

            query = """
                INSERT INTO node
                (name, display_name, primary_id_type, tax_id, type, pathways, source, function)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            """

            self.cursor.execute(query, (
                node_dict['name'],
                node_dict.get('display_name', ''),
                node_dict.get('primary_id_type'),
                node_dict.get('tax_id'),  # Can be None
                node_dict.get('type', 'protein'),
                node_dict.get('pathways', ''),
                node_dict.get('source', ''),
                node_dict.get('function', '')
            ))
            self.db.commit()

            node_dict['id'] = self.cursor.lastrowid
            self.check_node_dict_identifiers(node_dict)

        elif ('id' not in node_dict) and existing_node:
            node_dict['id'] = existing_node['id']

    def insert_unique_node(self, node_dict):
        # Add defaults for new fields
        if 'display_name' not in node_dict:
            node_dict['display_name'] = node_dict.get('name', '')

        if 'type' not in node_dict:
            node_dict['type'] = 'protein'  # Default type

        query = """
            INSERT INTO node
            (name, display_name, primary_id_type, tax_id, type, pathways, source, function)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """

        self.cursor.execute(query, (
            node_dict['name'],
            node_dict.get('display_name', ''),
            node_dict.get('primary_id_type'),
            node_dict.get('tax_id'),
            node_dict.get('type', 'protein'),
            node_dict.get('pathways', ''),
            node_dict.get('source', ''),
            node_dict.get('function', '')
        ))
        self.db.commit()

        node_dict['id'] = self.cursor.lastrowid
        self.check_node_dict_identifiers(node_dict)

    def get_node(self, node_name, node_tax_id=None):
        if node_tax_id is None:
            query = "SELECT * FROM node WHERE name = ?"
            tup = (node_name,)
        else:
            query = "SELECT * FROM node WHERE name = ? AND tax_id = ?"
            tup = (node_name, node_tax_id)

        self.cursor.execute(query, tup)
        self.db.commit()

        answer = self.cursor.fetchone()

        if answer:
            # Adjust to the new schema with display_name and type fields
            node_dict = {
                'id': answer[0],
                'name': answer[1],
                'primary_id_type': answer[2],
                'display_name': answer[3],
                'tax_id': answer[4],
                'type': answer[5],
                'pathways': answer[6],
                'source': answer[7],
                'function': answer[8]
            }

            # Also fetch identifiers
            node_dict.update(self.get_node_identifiers(node_dict['id']))

            return node_dict
        else:
            return None

    def get_node_by_name(self, node_name):
        query = "SELECT * FROM node WHERE name = ?"
        self.cursor.execute(query, (node_name,))
        self.db.commit()

        answer = self.cursor.fetchone()

        if answer:
            # Adjust to the new schema with display_name and type fields
            node_dict = {
                'id': answer[0],
                'name': answer[1],
                'primary_id_type': answer[2],
                'display_name': answer[3],
                'tax_id': answer[4],
                'type': answer[5],
                'pathways': answer[6],
                'source': answer[7],
                'function': answer[8]
            }

            # Also fetch identifiers
            node_dict.update(self.get_node_identifiers(node_dict['id']))

            return node_dict
        else:
            return None

    def get_node_by_id(self, id):
        self.cursor.execute("SELECT * FROM node WHERE id = ? ", (id, ))
        answer = self.cursor.fetchone()
        if not answer:
            return None
        else:
            # Adjust to the new schema with display_name and type fields
            id, name, primary_id_type, display_name, tax_id, mol_type, pathways, source, function = answer
            node_dict = {
                'id': answer[0],
                'name': answer[1],
                'primary_id_type': answer[2],
                'display_name': answer[3],
                'tax_id': answer[4],
                'type': answer[5],
                'pathways': answer[6],
                'source': answer[7],
                'function': answer[8]
            }
            # Also fetch identifiers
            node_dict.update(self.get_node_identifiers(id))

            return node_dict

    def update_node(self, node_dict):
        for k, v in node_dict.items():
            if not v and k != 'tax_id':  # Allow tax_id to be None for non-biological entities
                node_dict[k] = '-'

        # Handle identifier updates separately
        for id_type in ['kegg_id', 'uniprot_id', 'pubchem_id', "pubmed_id", "hgnc_id", "ensg_id"]:
            if id_type in node_dict:
                self.update_node_identifier(node_dict['id'], id_type, node_dict[id_type])
                del node_dict[id_type]  # Remove from main dict after handling

        tup = (
            node_dict['name'],
            node_dict.get('primary_id_type'),
            node_dict.get('display_name', node_dict['name']),
            node_dict.get('tax_id'),
            node_dict.get('type', 'protein'),
            node_dict.get('pathways', '-'),
            node_dict.get('source', '-'),
            node_dict.get('function', '-'),
            node_dict['id']
        )

        query = """
            UPDATE node
            SET name = ?, primary_id_type = ?, display_name = ?, tax_id = ?, type = ?,
                pathways = ?, source = ?, function = ?
            WHERE id = ?;
        """

        self.cursor.execute(query, tup)
        self.db.commit()

    def insert_edge(self, interactor_a_dict, interactor_b_dict, edge_dict):
        query = """
                INSERT INTO `edge` (
                `interactor_a_node_id`,
                `interactor_b_node_id`,
                `interactor_a_node_name`,
                `interactor_b_node_name`,
                `layer`,
                `source_db`,
                `interaction_types`
                )
                VALUES ( ?, ?, ?, ?, ?, ?, ?)
                """

        tup = (
            interactor_a_dict['id'],
            interactor_b_dict['id'],
            interactor_a_dict['name'],
            interactor_b_dict['name'],
            edge_dict['layer'],
            edge_dict['source_db'],
            edge_dict['interaction_types']
        )

        self.db.execute(query, tup)
        self.db.commit()

    def insert_node_identifier(self, node_id, id_type, id_value, is_primary=False):
        query = "INSERT INTO node_identifier (node_id, id_type, id_value, is_primary) VALUES (?, ?, ?, ?)"
        try:
            self.cursor.execute(query, (node_id, id_type, id_value, is_primary))
            self.db.commit()
        except sqlite3.IntegrityError:
            self.update_node_identifier(node_id, id_type, id_value, is_primary)

    def update_node_identifier(self, node_id, id_type, id_value, is_primary=False):
        self.cursor.execute("DELETE FROM node_identifier WHERE node_id = ? AND id_type = ?",
                            (node_id, id_type))
        # Insert new value if not empty
        if id_value and id_value != '-':
            self.cursor.execute("INSERT INTO node_identifier (node_id, id_type, id_value, is_primary) VALUES (?, ?, ?, ?)",
                               (node_id, id_type, id_value, is_primary))
        self.db.commit()

    def get_node_identifiers(self, node_id):
        """Get all identifiers for a node"""
        self.cursor.execute("SELECT id_type, id_value FROM node_identifier WHERE node_id = ?", (node_id,))
        identifiers = {}
        for id_type, id_value in self.cursor.fetchall():
            identifiers[id_type] = id_value
        return identifiers

    def get_node_by_identifier(self, id_type, id_value):
        """Find a node by one of its identifiers"""
        query = """
            SELECT n.* FROM node n
            JOIN node_identifier ni ON n.id = ni.node_id
            WHERE ni.id_type = ? AND ni.id_value = ?
        """
        self.cursor.execute(query, (id_type, id_value))
        answer = self.cursor.fetchone()

        if not answer:
            return None

        # Adjust to the new schema with display_name and type fields
        id, name, primary_id_type, display_name, tax_id, mol_type, pathways, source, function = answer
        node_dict = {
            "id": id,
            "name": name,
            "primary_id_type": primary_id_type,
            "display_name": display_name,
            "tax_id": tax_id,
            "type": mol_type,
            "pathways": pathways,
            "source": source,
            "function": function
        }

        # Also fetch all identifiers
        node_dict.update(self.get_node_identifiers(id))

        return node_dict

    def save_db_to_file(self, db_file_name):
        if '.db' not in db_file_name:
            export_file = db_file_name + '.db'
        else:
            export_file = db_file_name

        self.create_db(export_file)

        db_file_name = db_file_name.split('/')[-1]
        db_name = db_file_name.replace(".db", "")

        tup = (export_file, db_name)
        self.db.execute("ATTACH ? as ?", tup)
        self.db.commit()

        # Copy node table
        attached_table = (db_name + '.node')
        self.db.execute("INSERT INTO %s SELECT * FROM node" % attached_table)
        self.db.commit()

        # Copy edge table
        attached_table = (db_name + '.edge')
        self.db.execute("INSERT INTO %s SELECT * FROM edge" % attached_table)
        self.db.commit()

        # Copy node_identifier table
        try:
            attached_table = (db_name + '.node_identifier')
            self.db.execute("INSERT INTO %s SELECT * FROM node_identifier" % attached_table)
            self.db.commit()
        except sqlite3.OperationalError:
            # Table might not exist in the target database
            pass

        self.db.close()
