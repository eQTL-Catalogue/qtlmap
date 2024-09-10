#!/usr/bin/env python

import duckdb
import argparse

argparser = argparse.ArgumentParser()
argparser.add_argument('-v', help='Variant info file from VCF.', required=True)
argparser.add_argument('-o', help='The output pq file name', required=True)

args = argparser.parse_args()

variant_info_from_vcf = args.v
output = args.o

con = duckdb.connect()

con.execute(f"""
    COPY (
        SELECT
            chromosome, position, variant, ref, alt, type, ac, an, r2
        FROM read_csv_auto('{variant_info_from_vcf}', delim='\t', compression='gzip', header=False, 
            columns={{'chromosome': 'VARCHAR', 'position': 'INTEGER', 'variant': 'VARCHAR', 'ref': 'VARCHAR', 
                            'alt': 'VARCHAR', 'type': 'VARCHAR', 'ac': 'INTEGER', 'an': 'INTEGER', 'maf': 'DOUBLE','r2': 'VARCHAR'}})) 
    TO '{output}' (FORMAT 'parquet');
""")
con.close()