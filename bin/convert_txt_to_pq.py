#!/usr/bin/env python

import duckdb
import argparse
import json 

argparser = argparse.ArgumentParser()
argparser.add_argument('-i', help='Input file (optional,None if tmp is missing)', required=False, default=None)
argparser.add_argument('-t', help='tmp', required=False, default=0, type=int)
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

if input_file and not args.tmp:
    con.execute(f"""
        COPY (
            SELECT
                {columns}
            FROM read_csv_auto('{input_file}'{schema_query})) 
        TO '{output}' (FORMAT 'parquet');
    """)
elif not input_file and args.tmp:
    con.execute(f"""
        CREATE TABLE temp_table (
            phenotype_id VARCHAR,
            median_tpm DOUBLE 
        );
        COPY temp_table TO '{args.o}' (FORMAT 'parquet');
    """)

con.close()