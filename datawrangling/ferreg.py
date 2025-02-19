import pandas as pd
from pathlib import Path
import sqlite3


def parse_table_name(s: str) -> str:
    ss = s.split("/")
    return ss[-2]+"|"+ss[-1].split(".")[0]


ferreg_dir = Path("./sources/ferreg/")
ferreg_tables = list(ferreg_dir.glob('**/*.*'))

cols = []
with sqlite3.connect('ferreg.db') as conn:
    for table in ferreg_tables:
        try:
            k = parse_table_name(str(table))
            df = pd.read_csv(table, delimiter='\t')
            v = df.columns.tolist
            cols.append({k: v})
            df.to_sql(k, conn, if_exists='replace', index=False)
        except Exception as e:
            print(f"Something went wrong with {table}: {e}")
            continue
