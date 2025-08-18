from config import OUTPUTS_DIR
from database.external_db import DBconnector
import matplotlib.pyplot as plt

ferreg = DBconnector(OUTPUTS_DIR / "ferreg.db")
giga_query = """
SELECT
    hub.unique_id,
    COALESCE(t.uniprot_id, hub.target_id) AS target_id,
    COALESCE(r.External_id, hub.regulator_id) AS regulator_id,
    COALESCE(d."Disease ICD", hub.disease_id) AS disease_id,
    COALESCE(l.drug_name, hub.drug_id) AS drug_id,
    reg.*
FROM target_regulator_drug_disease_pair AS hub
LEFT JOIN regulation_information AS reg
    ON hub.unique_id = reg.unique_id
LEFT JOIN general_regulator AS r
    ON hub.regulator_id = r.regulator_id
LEFT JOIN general_target AS t
    ON hub.target_id = t.target_id
LEFT JOIN general_disease AS d
    ON hub.disease_id = d.disease_id
LEFT JOIN general_drug AS l
    ON hub.drug_id = l.drug_id
"""
res = ferreg.query_to_dataframe(giga_query)
res.to_csv("./giga_query_res.csv")

smart_query = """
SELECT
    disease_id,
    GROUP_CONCAT(regulator_id, ',') as regulators,
    COUNT(regulator_id) as reg_count
FROM target_regulator_drug_disease_pair
WHERE regulator_id != "."
GROUP BY disease_id
"""
grouped = ferreg.query_to_dataframe(smart_query)
grouped.sort_values("reg_count", ascending=False)
grouped.reg_count.unique()
grouped.plot.bar("disease_id", "reg_count")
plt.show()
