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

    def insert_node(self, node_dict):

        if ('id' not in node_dict) and (not self.get_node(node_dict['name'], node_dict['tax_id'])):
            query = "INSERT INTO node (name, gene_name, tax_id, pathways) VALUES (?, ?, ?, ?)"

            self.cursor.execute(query, (
                node_dict['name'],
                node_dict['gene_name'],
                node_dict['tax_id'],
                node_dict['pathways']
            ))
            self.db.commit()

            node_dict['id'] = self.cursor.lastrowid

        elif ('id' not in node_dict) and (self.get_node(node_dict['name'], node_dict['tax_id'])):
            node_dict['id'] = self.get_node(node_dict['name'], node_dict['tax_id'])['id']

    def insert_unique_node(self, node_dict):

        query = "INSERT INTO node (name, gene_name, tax_id, pathways) VALUES (?, ?, ?, ?)"

        self.cursor.execute(query, (
            node_dict['name'],
            node_dict['gene_name'],
            node_dict['tax_id'],
            node_dict['pathways']))
        self.db.commit()

        node_dict['id'] = self.cursor.lastrowid

    def get_node(self, node_name, node_tax_id):
        query = "SELECT * FROM node WHERE name = ? AND tax_id = ?"
        tup = (node_name, node_tax_id)

        self.cursor.execute(query, tup)
        self.db.commit()

        answer = self.cursor.fetchone()

        if answer:
            node_dict = {
                'id': answer[0],
                'name': answer[1], #primary id of the database
                'gene_name': answer[2], #other ids in mi: format
                'tax_id': answer[3],
                'pathways': answer[4], #only in 0-1st layers
            }
            return node_dict
        else:
            return None

    def get_node_by_id(self,id):
        self.cursor.execute("SELECT * FROM node WHERE id = ? ", (id, ))
        answer = self.cursor.fetchone()
        if not answer:
            return None
        else:
            if answer:
                id, name, gene_name, tax_id, pathways = answer
                node_dict = {
                    "id": id,
                    "name": name,
                    "gene_name": gene_name,
                    "tax_id": tax_id,
                    "pathways": pathways
                }
            else:
                return None
            return node_dict

    def update_node(self, node_dict):
        for k, v in node_dict.items():
            if not v:
                node_dict[k] = '-'

        tup = (node_dict['name'],
               node_dict['gene_name'],
               node_dict['tax_id'],
               node_dict['pathways'],
               node_dict['id'])

        query = """
            UPDATE node
            SET name = ?, gene_name = ?, tax_id = ?, pathways = ?
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
                `source_db`,
                `interaction_identifiers`,
                `interaction_types`
                )
                VALUES ( ?, ?, ?, ?, ?, ?, ?)
                """

        tup = (
            interactor_a_dict['id'],
            interactor_b_dict['id'],
            interactor_a_dict['name'],
            interactor_b_dict['name'],
            edge_dict['source_db'],
            edge_dict['interaction_identifiers'],
            edge_dict['interaction_types']
        )

        self.db.execute(query, tup)
        self.db.commit()

    def save_db_to_file(self, db_file_name):

        if '.db' not in db_file_name:
            export_file = db_file_name + '.db'
        else:
            export_file = db_file_name

#        file_db = self.create_db(export_file)
        self.create_db(export_file)

        db_file_name = db_file_name.split('/')[-1]

        db_name = db_file_name.replace(".db", "")

        tup = (export_file, db_name)
        self.db.execute("ATTACH ? as ?", tup)
        self.db.commit()

        attached_table = (db_name + '.node')
        self.db.execute("INSERT INTO %s SELECT * FROM node" % attached_table)
        self.db.commit()

        attached_table = (db_name + '.edge')
        self.db.execute("INSERT INTO %s SELECT * FROM edge" % attached_table)
        self.db.commit()

        self.db.close()
