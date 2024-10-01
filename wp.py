import pandas as pd
from Bio.KEGG.KGML.KGML_parser import read
import pickle
import apicalls.api_calls as api

pathway = read(open("./sources/kegg/hsa04216.xml", "r"))

kegg_gene_dict = dict()
for k, v in pathway.entries.items():
    kegg_gene_dict[k] = [v.name]
df_kegg = pd.DataFrame.from_dict(kegg_gene_dict)

all_ids = df_kegg.iloc[0].to_list()
ids = []
for i in all_ids:
    if 'hsa' in i and 'path' not in i:
        ids.append(i)
ids
db_dict = {'uniprot': "UniProtKB-Swiss-Prot",
           'kegg': "KEGG",
           'genewiki': "GeneWiki",
           'ensembl': "Ensembl",
           'gene': "Gene_Name"}
wp_id_dict = dict()
for i in ids:
    job_id = api.submit_uni_id_mapping(
        fromDB=db_dict["kegg"], toDB=db_dict['uniprot'], ids=[i], human=False
    )
    results = api.get_uni_id_mapping_results(job_id)
    print(results['results'][0]['to']['primaryAccession'])
    wp_id_dict[i] = results['results'][0]['to']['primaryAccession']
pickle.dump(wp_id_dict, open('kegg_id.p', 'wb'))
wp_id_dict
df_go = pd.read_csv("./sources/go/hsa_gene.csv", delimiter=',')
go_dict = dict()
for idr, row in df_go.iterrows():
    go_dict[row.bioentity.split(':')[1]] = row.bioentity_label
pickle.dump(go_dict, open('go_id.p', 'wb'))
