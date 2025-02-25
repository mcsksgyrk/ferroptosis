import pandas as pd
from config import OUTPUTS_DIR, SOURCES_DIR
import sqlite3


def make_sql_db(name: str):
    source_dir = SOURCES_DIR / name
    source_tables = list(source_dir.glob('**/*.*'))
    output_name = name+".db"
    cols = []
    with sqlite3.connect(OUTPUTS_DIR / output_name) as conn:
        for table in source_tables:
            if ".DS_Store" in str(table):
                continue
            try:
                k = str(table).split("/")[-1]
                df = pd.read_csv(table, sep=None, engine="python")
                v = df.columns.tolist
                cols.append({k: v})
                df.to_sql(k, conn, if_exists='replace', index=False)
            except Exception as e:
                print(f"Something went wrong with {table}: {e}")
                continue


make_sql_db("ferrdb")
make_sql_db("ferreg")
