from config import OUTPUTS_DIR
from database.external_db import DBconnector


db_path = OUTPUTS_DIR / "ferr_test.db"
db = DBconnector(db_path)
