import pandas as pd
from config import SOURCES_DIR, OUTPUTS_DIR
import sqlite3
from typing import List, Union, Optional, Any


class DBconnector:
    def __init__(self, db_path: str = str(OUTPUTS_DIR / "ferreg.db")):
        self.db_path = db_path

    def _get_connection(self):
        return sqlite3.connect(self.db_path)

    def get_table_names(self):
        with self._get_connection() as conn:
            cursor = conn.cursor()
            res = cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables = res.fetchall()
        return [row[0] for row in tables]

    def get_columns(self, table_name):
        with self._get_connection() as conn:
            cursor = conn.cursor()
            res = cursor.execute(f"PRAGMA table_info({table_name})")
            columns = res.fetchall()
        return columns

    def query(self, what: str, table: str, where: Optional[str] = None,
              values: Optional[List] = None) -> List:
        with self._get_connection() as conn:
            if values is not None and where is not None:
                placeholders = ",".join("?" * len(values))
                query = f'SELECT {what} FROM {table} WHERE {where} IN ({placeholders})'
                cursor = conn.cursor()
                res = cursor.execute(query, values).fetchall()
            else:
                query = f'SELECT {what} FROM {table}'
                cursor = conn.cursor()
                res = cursor.execute(query).fetchall()

            return [row[0] for row in res]

    def custom_query(self, query: str) -> List:
        with self._get_connection() as conn:
            cursor = conn.cursor()
            res = cursor.execute(query).fetchall()
            if len(res[0]) == 1:
                return [row[0] for row in res]
            else:
                return res

    def query_to_dataframe(self, query: str) -> pd.DataFrame:
        with self._get_connection() as conn:
            return pd.read_sql_query(query, conn)

    @staticmethod
    def compile_query_string(columns: str, table: str,
                            filters: Optional[List[str]] = None) -> str:
        query = f"SELECT {columns} FROM {table}"
        if filters:
            query += " WHERE " + " AND ".join(filters)
        return query
