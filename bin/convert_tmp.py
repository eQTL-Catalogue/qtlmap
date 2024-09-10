#!/usr/bin/env python

import duckdb
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('-m', help='Median tmp file', required=True)
argparser.add_argument('-o', help='The output pq file name', required=True)

args = argparser.parse_args()

median_tmp_file = args.m
output = args.o

con = duckdb.connect()

con.execute(f"""
    COPY (
        SELECT
            phenotype_id, median_tpm
        FROM read_csv_auto('{median_tmp_file}', delim='\t', compression='gzip', header=True)) 
    TO '{output}' (FORMAT 'parquet');
""")
con.close()