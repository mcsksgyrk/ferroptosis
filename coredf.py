import pandas as pd
import pickle
import os
import apicalls.api_calls as api

kegg = pickle.load(open('kegg_id.p', 'rb'))
go = pickle.load(open('go_id.p', 'rb'))
wp = pickle.load(open('wp_id.p', 'rb'))

go_vals = set(go.values())
kegg_vals = set(kegg.values())
wp_vals = set(wp.values())
unique_prs = go_vals | kegg_vals | wp_vals
ids = list(unique_prs)

db_dict = {'uniprot_to': "UniProtKB-Swiss-Prot",
           'kegg': "KEGG",
           'genewiki': "GeneWiki",
           'ensembl': "Ensembl",
           'gene': "Gene_Name",
           "uniprot_from": "UniProtKB_AC-ID"}

res, err = api.convert_from_uniprot_id(db_dict["gene"], ids, human=False)

unip = res.keys()
gene = res.values()
to_csv = pd.DataFrame({'uniprotID': unip, 'gene': gene})
to_csv.to_csv('prs_from_kegg_go_wp.csv',index=False)
