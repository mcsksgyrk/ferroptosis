import pandas as pd
from parsers.source_parsers import KEGGPathwayParser, WikiPathwayParser
from apicalls.api_oop import UniProtClient
from typing import Dict, List
from pathlib import Path
from config import OUTPUTS_DIR


def make_tables(data: Dict[int, List[str]], source: str, core: int) -> pd.DataFrame:
    df = pd.DataFrame(list(data.items()), columns=["source_id", "uniprot_id"])
    df['source'] = source
    df['is_core'] = core
    return df


def make_go_table(path: Path, core: int) -> pd.DataFrame:
    go_df = pd.read_csv(path)
    standard_df = pd.DataFrame({
        "uniprot_id": go_df.bioentity.str.replace("UniProtKB:", ""),
        "source": "GO",
        "is_core": core
    })
    return standard_df


def combine_sources(dataframes_list: list[pd.DataFrame]) -> pd.DataFrame:
    combined_df = pd.concat(dataframes_list)
    grouped = combined_df.groupby(['uniprot_id', 'is_core'])
    result_df = grouped.agg({
        'source': lambda x: ';'.join(sorted(set(x)))
    }).reset_index()
    return result_df


kegg_parser = KEGGPathwayParser()
wikipath_parser = WikiPathwayParser()
wp_src = wikipath_parser.read_pathway("WP4313-datanodes.tsv")
kegg_src = kegg_parser.read_pathway("hsa04216.xml")
wp_id_list = wikipath_parser.extract_gene_ids(wp_src)
kegg_id_list = kegg_parser.extract_gene_ids(kegg_src)

uniprot_client = UniProtClient()
f, g = uniprot_client.convert_to_uniprot_id("KEGG", kegg_id_list, False)
print(f"{g} couldn't be mapped to uniprot id")
kegg_df = make_tables(f, "kegg", 1)
wp_df = make_tables(wp_id_list, "wikipathways", 1)
go_df = make_go_table("sources/go/hsa_gene.csv", 1)
result_df = combine_sources([wp_df, kegg_df, go_df])
result_df.to_csv(OUTPUTS_DIR / "core_geneProducts.csv", index=False)
