""" Module to set up database connection """

import mysql.connector

class database_connection():

    def __init__(self,user="af",password="afPass",host="localhost",port="3308",database="atlas_moving_objects"):
        self.config = {
          'user': user,
          'password': password,
          'host': host,
          'port': port,
          'database': database,
          'raise_on_warnings': True
        }

    def connect(self):
        cnx = mysql.connector.connect(**self.config)
        return cnx

    def disconnect(self):
        self.cnx.close()

# if __name__ == "__main__":
#     cnx=database_connection().connect()
#     print(cnx)
#     cnx.disconnect()
