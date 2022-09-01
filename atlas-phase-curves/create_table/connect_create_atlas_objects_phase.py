"""
"""

import mysql.connector
import pandas as pd
from optparse import OptionParser

# connect to the database
config = {
  'user': 'af',
  'password': 'afPass',
  'host': 'localhost',
  'port': '3308',
  'database': 'atlas_moving_objects',
  # 'database': 'jamie_test_db',
  'raise_on_warnings': True
}
cnx = mysql.connector.connect(**config)
cursor =cnx.cursor()

# sql_file = "create_atlas_phase_fits.sql"
sql_file = "create_atlas_phase_fits_app.sql"

# read in the sql query (note that this is one query split over multiple lines)
# multiple queries will probably require each line to be executed separately
with open(sql_file) as f:
    qry = f.readlines()
qry="".join(qry)
print(qry)

# perform the query
cursor.execute(qry)
cnx.commit()

cnx.close()
