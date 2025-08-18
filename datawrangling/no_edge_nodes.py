from config import OUTPUTS_DIR, PROJECT_ROOT
from database.external_db import DBconnector

fer_path = OUTPUTS_DIR / 'test_omnipath.db'
ferrdb = DBconnector(fer_path)

lonely_query = """
SELECT n.id, n.name, n.source
FROM node AS n
WHERE n.type = 'protein'
AND n.id NOT IN (
    SELECT DISTINCT interactor_a_node_id FROM edge
    UNION
    SELECT DISTINCT interactor_b_node_id FROM edge
)
"""
res = ferrdb.custom_query(lonely_query)
len(res)
