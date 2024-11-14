import pandas as pd
import numpy as np
import apicalls.api_calls as api

df_driver = pd.read_csv('./sources/ferrdb/ferroptosis_early_preview_upto20221231/driver.csv')
df_supp = pd.read_csv('./sources/ferrdb/ferroptosis_early_preview_upto20221231/suppressor.csv')
df_mark = pd.read_csv('./sources/ferrdb/ferroptosis_early_preview_upto20221231/marker.csv')
column_names = df_supp.columns

driver_gene_list = df_driver[(df_driver.Confidence=="Validated") & (df_driver.Exp_organism.str.contains("Human"))].UniProtAC.unique()
suppr_gene_list = df_supp[(df_supp.Confidence=="Validated") & (df_supp.Exp_organism.str.contains("Human"))].UniProtAC.unique()
mark_gene_list = df_mark[(df_mark.Confidence=="Validated") & (df_mark.Exp_organism.str.contains("Human"))].UniProtAC.unique()

interact_genes = pd.read_csv("prs_from_kegg_go_wp.csv")
connceter_prs = interact_genes.uniprotID.to_list()

len(connceter_prs)

allin1 = np.concatenate((suppr_gene_list, driver_gene_list, mark_gene_list))
unique_ferrdb_prs = np.unique(allin1.astype(str))
missing = []
for pr in connceter_prs:
    if pr not in unique_ferrdb_prs:
        missing.append(pr)

res1, rip1 = api.convert_from_uniprot_id("Gene_Name", missing, human=False)
res2, rip2 = api.convert_from_uniprot_id("Gene_Name", allin1, human=False)

all_prs = res2 | res1
core_prs_df = pd.DataFrame(all_prs.items(), columns=["UniProtAC", "Gene_Name"])
kek = core_prs_df.UniProtAC.to_list()
missing_markers = []
for pr in mark_gene_list:
    if pr not in kek:
        missing_markers.append(pr)
res3, rip3 = api.convert_from_uniprot_id("Gene_Name", missing_markers, human=False)
all_prs = res3 | res2 | res1
core_prs_df = pd.DataFrame(all_prs.items(), columns=["UniProtAC", "Gene_Name"])
core_prs_df.to_csv("core_proteins.csv")
core_prs_df
