import sqlite3
import pandas as pd
from config import OUTPUTS_DIR
from typing import Dict, Set
import json



merged_db_path = OUTPUTS_DIR / 'final_network_with_omnipath.db'

conn = sqlite3.connect(merged_db_path)

fer_layer1_query = """
    SELECT e.interactor_a_node_name, e.interactor_b_node_name, e.layer
    FROM edge e
    JOIN node n1 ON e.interactor_a_node_name = n1.name
    JOIN node n2 ON e.interactor_b_node_name = n2.name
    WHERE (e.layer = "0.0" OR e.layer = "0")
    AND e.source_db NOT LIKE '%ARN%'
    AND n1.primary_id_type = "uniprot_id"
    AND n2.primary_id_type = "uniprot_id"
"""
layer_1_df = pd.read_sql_query(fer_layer1_query, conn)
fer_layer1 = set(layer_1_df.interactor_a_node_name)
fer_layer0 = set(layer_1_df.interactor_b_node_name)
len(fer_layer0)
