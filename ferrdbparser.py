import pandas as pd

df_driver = pd.read_csv("./sources/ferrdb/ferroptosis_early_preview_upto20221231/driver.csv")

df_driver.Confidence.unique()
