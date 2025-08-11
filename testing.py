import pandas as pd
from pathlib import Path
from config import OUTPUTS_DIR, SOURCES_DIR, PROJECT_ROOT
import sqlite3

scrna = SOURCES_DIR / "cellmaker"
dbpath = OUTPUTS_DIR / "test_omnipath.db"

df = pd.read_excel(scrna / "Cell_marker_Human.xlsx")
df.columns
conn = sqlite3.connect(dbpath)
nodes = pd.read_sql_query("SELECT * FROM node WHERE type = 'protein'", conn)
node_ids = nodes['name'].tolist()
final_df = df[df['UNIPROTID'].isin(node_ids)]
len(final_df.UNIPROTID.unique())
final_df.cellontology_id.unique()
final_df.uberonongology_id.unique()
