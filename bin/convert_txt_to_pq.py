#!/usr/bin/env python

import duckdb
import argparse
import json 

argparser = argparse.ArgumentParser()
argparser.add_argument('-i', help='Input file', required=True)
argparser.add_argument('-c', '--columns', required=True, help="List of file columns.")
argparser.add_argument('-o', help='The output pq file name', required=True)
argparser.add_argument('-s', '--schema', help='Query ending as json dict of column types',type=str,default="", required=False)

args = argparser.parse_args()

input_file = args.i
output = args.o
columns = args.columns
schema_query = ""
if args.schema:
    columns_dict = json.loads(args.schema)
    schema_query = f", columns={columns_dict}"
con = duckdb.connect()

con.execute(f"""
    COPY (
        SELECT
            {columns}
        FROM read_csv_auto('{input_file}'{schema_query})) 
    TO '{output}' (FORMAT 'parquet');
""")
con.close()