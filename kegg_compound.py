from config import SOURCES_DIR


compounds = set()
f_path = SOURCES_DIR / "kegg/kegg_compounds.txt"
with open(f_path, 'r') as f:
    for line in f:
        kegg_entry = line.strip().split('\t', 1)[1]
        names = [name.strip().lower() for name in kegg_entry.split(';')]
        compounds.update(names)
