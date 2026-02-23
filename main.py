from datawrangling.make_sql_dbs import sources_to_sql_schema
from datawrangling.transform_core import convert_kegg_source
from datawrangling.transform_ferrdb import convert_ferrdb_source
from datawrangling.transform_ferreg import convert_ferreg_source
from database.merger import merger_sources
from database.merger_disease import migrate_metadata
from datawrangling.edges_from_omnipath import extend_merged_db_with_omnipath

# converts the source files to SQL schemas
sources_to_sql_schema()
# convert sql data into unified SQL schema
convert_kegg_source()
convert_ferrdb_source()
convert_ferreg_source()
# merges the unifiely formatted datafiles
merger_sources()
migrate_metadata()
# looks for new edges in OmniPath database
extend_merged_db_with_omnipath()
