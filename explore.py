import pandas as pd
from Bio.KEGG.KGML.KGML_parser import read
import pickle
import apicalls.api_calls as api

pathway = read(open("./sources/kegg/hsa04216.xml", "r"))
pathway.entries
df_wp = pd.read_csv('./sources/wikipathways/WP4313-datanodes.tsv',
                    delimiter='\t')
genes_wp = df_wp[df_wp.Type == "GeneProduct"]
ids = genes_wp.Label.to_list()
ids
db_dict = {'uniprot': "UniProtKB-Swiss-Prot",
           'kegg': "KEGG",
           'genewiki': "GeneWiki",
           'ensembl': "Ensembl",
           'gene': "Gene_Name"}
wp_id_dict = dict()
for i in ids:
    job_id = api.submit_uni_id_mapping(
        fromDB=db_dict["gene"], toDB=db_dict["uniprot"], ids=[i]
    )
    results = api.get_uni_id_mapping_results(job_id)
    print(results['results'][0]['to']['primaryAccession'])
    wp_id_dict[i] = results['results'][0]['to']['primaryAccession']

pickle.dump(wp_id_dict, open('wp_id.p', 'wb'))

pathway = read(open("./sources/kegg/hsa04216.xml", "r"))

kegg_gene_dict = dict()
for k, v in pathway.entries.items():
    kegg_gene_dict[k] = [v.name]
df_kegg = pd.DataFrame.from_dict(kegg_gene_dict)

all_ids = df_kegg.iloc[0].to_list()
ids = []
for i in all_ids:
    if 'hsa' in i and 'path' not in i:
        id_ = i.split(' ')
        for j in id_:
            ids.append(j)
ids
db_dict = {'uniprot': "UniProtKB-Swiss-Prot",
           'kegg': "KEGG",
           'genewiki': "GeneWiki",
           'ensembl': "Ensembl",
           'gene': "Gene_Name",
           'uniprot_from': "UniProtKB_AC-ID"}
rippers = []
for i in ids:
    job_id = api.submit_uni_id_mapping(
        fromDB=db_dict["uniprot_from"], toDB=db_dict["kegg"], ids=[i], human=False
    )
    results = api.get_uni_id_mapping_results(job_id)
    print(results)
    try:
    except IndexError:
        print(f'ripped id {i} with job id {job_id}')
        rippers.append(i)

pickle.dump(wp_id_dict, open('kegg_id.p', 'wb'))
len(wp_id_dict)
df_go = pd.read_csv("./sources/go/hsa_gene.csv", delimiter=',')
go_dict = dict()
for idr, row in df_go.iterrows():
    go_dict[row.bioentity_label] = row.bioentity.split(':')[1]
pickle.dump(go_dict, open('go_id.p', 'wb'))
