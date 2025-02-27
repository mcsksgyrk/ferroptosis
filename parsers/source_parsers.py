from pathlib import Path
from typing import Dict, List
import logging
import pandas as pd
from config import SOURCES_DIR, OUTPUTS_DIR


class PathwayParser:
    def __init__(self, data_dir: Path, logger_name=None):
        self.data_dir = Path(data_dir)
        logger_name = logger_name or __name__
        self.logger = logging.getLogger(logger_name)


class KEGGPathwayParser(PathwayParser):
    def __init__(self, data_dir: Path = SOURCES_DIR / "kegg"):
        super().__init__(data_dir=data_dir)

    def _find_names(self, line: str) -> List[str]:
        name_start = line.find('name="') + 6
        name_end = line.find('"', name_start)
        name_value = line[name_start:name_end]
        return name_value.split()

    def _find_id(self, line: str) -> int:
        id_start = line.find('id="') + 4
        id_end = line.find('"', id_start)
        return int(line[id_start:id_end])

    def read_pathway(self, filename: str) -> Dict[int, List[str]]:
        try:
            entries = {}
            with open(self.data_dir / filename) as f:
                for line in f:
                    if "<entry" in line:
                        k = self._find_id(line)
                        v = self._find_names(line)
                        entries[k] = v
            return entries
        except FileNotFoundError:
            self.logger.error(f"Pathway file not found: {filename}")
            raise
        except Exception as e:
            self.logger.error(f"Error reading pathway file: {str(e)}")
            raise

    def extract_gene_ids(self, pathway_dict: Dict[int, List[str]]) -> List[str]:
        try:
            human_genes = []
            for _, values in pathway_dict.items():
                for item in values:
                    if item.startswith("hsa") and "path" not in item:
                        human_genes.append(item)
            self.logger.info(f"Extracted {len(human_genes)} gene IDs")
            return human_genes
        except Exception as e:
            self.logger.error(f"Error extracting gene IDs: {str(e)}")
            raise


class WikiPathwayParser(PathwayParser):
    def __init__(self, data_dir: Path = SOURCES_DIR / "wikipathways"):
        super().__init__(data_dir=data_dir)

    def read_pathway(self, filename: str) -> pd.DataFrame:
        try:
            filepath = self.data_dir / filename
            df = pd.read_csv(filepath, delimiter='\t')
            self.logger.info(f"Read {len(df)} data nodes from {filename}")
            return df
        except Exception as e:
            self.logger.error(f"Error reading pathway file: {str(e)}")
            raise

    def _get_canonical_id(self, id_list: List[str]) -> str:
        cans = ['O', 'Q', 'P']
        try:
            for id in id_list:
                if id[0] in cans:
                    return id
        except Exception as e:
            self.logger.info(f"Error {e}")
            raise

    def extract_gene_ids(self, df: pd.DataFrame, canonical = True) -> Dict[int, list[str]]:
        try:
            uniprot_ids = {}
            geneProducts = df[df.Type == "GeneProduct"]
            for idx, row in geneProducts.iterrows():
                k = int(idx)
                if canonical:
                    ax = row.UniProt.replace("uniprot:", "").split(";")
                    v = self._get_canonical_id(ax)
                else:
                    v = row.UniProt.replace("uniprot:", "").split(";")
                uniprot_ids[k] = v
            self.logger.info(f"Extracted {len(uniprot_ids)} gene IDs")
            return uniprot_ids
        except Exception as e:
            self.logger.error(f"Error: {str(e)}")
            raise
