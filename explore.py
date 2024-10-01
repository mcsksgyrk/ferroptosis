import pandas as pd
from Bio.KEGG.KGML.KGML_parser import read
import pickle
import apicalls.api_calls as api
import json

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
