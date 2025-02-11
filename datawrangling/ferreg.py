import pandas as pd
from typing import List, Dict, Union, TypedDict
from pathlib import Path
import sqlite3
import os


def parse_table_name(s: str) -> str:
    ss = s.split("/")
    return ss[-2]+"|"+ss[-1].split(".")[0]


ferreg_dir = Path("./sources/ferreg/")
ferreg_tables = list(ferreg_dir.glob('**/*.*'))

conn = sqlite3.connect('ferreg.db')
cols = []
for table in ferreg_tables:
    k = parse_table_name(str(table))
    df = pd.read_csv(table, delimiter='\t')
    v = df.columns.tolist
    cols.append({k: v})
    df.to_sql(k, conn, if_exists='replace', index=False)

conn.close()
