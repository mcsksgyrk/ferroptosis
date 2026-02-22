from datawrangling.make_sql_dbs import sources_to_sql_schema
from datawrangling.transform_core import convert_kegg_source
from datawrangling.transform_ferrdb import convert_ferrdb_source
from datawrangling.transform_ferreg import convert_ferreg_source
from database.merger import merger_sources


sources_to_sql_schema()
convert_kegg_source()
convert_ferreg_source()
convert_ferrdb_source()
merger_sources()
