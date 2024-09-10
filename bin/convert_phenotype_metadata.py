#!/usr/bin/env python

import duckdb
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('-m', help='Phenotype metadata file', required=True)
argparser.add_argument('-o', help='The output pq file name', required=True)

args = argparser.parse_args()

phenotype_metadata = args.m
output = args.o

con = duckdb.connect()

con.execute(f"""
    COPY (
        SELECT
            phenotype_id, group_id, gene_id
        FROM read_csv_auto('{phenotype_metadata}', delim='\t', header=True)) 
    TO '{output}' (FORMAT 'parquet');
""")
con.close()