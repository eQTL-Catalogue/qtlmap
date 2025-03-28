#!/usr/bin/env python

import argparse
import duckdb
import time


argparser = argparse.ArgumentParser()
argparser.add_argument('-v', help='The variant info file, extracted from vcf', required=True)
argparser.add_argument('-r', help='The variant_id rsid map file', required=True)
argparser.add_argument('-s', help='The summary_statistics file (nominal run output)', required=True)
argparser.add_argument('-p', help='The phenotype metadata file', required=True)
argparser.add_argument('-m', help='The median_tpm file', required=False, default=None)
argparser.add_argument('-o', help='The output file name', required=True)
argparser.add_argument('-a', help='Starting position', required=True,type=int)
argparser.add_argument('-b', help='Ending position', required=True,type=int)
argparser.add_argument('-t', help='tpm_missing', required=True,type=int)
argparser.add_argument('-e',required=True, type=str, help="Memory limit in GB for DuckDB.")




args = argparser.parse_args()

rsid_map_file = args.r
variant_info_from_vcf = args.v
nominal_run_output = args.s
phenotype_metadata = args.p
median_tmp_file = args.m
output_file = args.o
start_pos = args.a
end_pos = args.b
tpm_file_missing = args.t
memory_limit = args.e




def main(memory_limit):
    memory_limit = f"{float(memory_limit) * 0.8:.1f}"
    start = time.time()
    con = duckdb.connect()
    if tpm_file_missing:
        median_tpm_query = ""
        median_tpm_column = "NULL AS median_tpm" 
    else:
        median_tpm_query = f"""
            LEFT JOIN (
                SELECT phenotype_id, median_tpm 
                FROM read_parquet('{median_tmp_file}')
            ) AS mt ON nro.molecular_trait_id = mt.phenotype_id
        """
        median_tpm_column = "mt.median_tpm"
    con.execute(f"""
        SET memory_limit='{memory_limit}GB';
        COPY (
            SELECT 
                nro.molecular_trait_id, 
                vi.chromosome, 
                vi.position, 
                vi.ref, 
                vi.alt, 
                nro.variant, 
                nro.ma_samples, 
                nro.maf, 
                nro.pvalue, 
                nro.beta, 
                nro.se, 
                vi.type, 
                vi.ac, 
                vi.an, 
                vi.r2, 
                pm.group_id AS molecular_trait_object_id, 
                pm.gene_id, 
                {median_tpm_column}, 
                rm.rsid
            FROM read_csv_auto('{nominal_run_output}', delim='\t', header=False, 
                    columns={{'molecular_trait_id': 'VARCHAR', 'variant': 'VARCHAR', 'X': 'VARCHAR', 
                            'ma_samples': 'INTEGER', 'ma_count': 'INTEGER', 'maf': 'DOUBLE', 
                            'pvalue': 'DOUBLE', 'beta': 'DOUBLE', 'se': 'DOUBLE'}}) AS nro
            -- INNER JOIN on variant info table (only keep rows where variants match)
            INNER JOIN (
                SELECT 
                    chromosome, 
                    position, 
                    variant, 
                    ref, 
                    alt, 
                    type, 
                    ac, 
                    an, 
                    r2
                FROM read_parquet('{variant_info_from_vcf}')
            ) AS vi ON nro.variant = vi.variant
            LEFT JOIN (
                SELECT 
                    variant, 
                    rsid, 
                    position
                FROM read_parquet('{rsid_map_file}')
                WHERE position BETWEEN {start_pos} AND {end_pos}
            ) AS rm ON nro.variant = rm.variant
            {median_tpm_query}
            LEFT JOIN (
                SELECT phenotype_id, group_id, gene_id 
                FROM read_parquet('{phenotype_metadata}')
            ) AS pm ON nro.molecular_trait_id = pm.phenotype_id)
                TO '{output_file}' (FORMAT 'parquet');
        """)      
    con.close()              
    end = time.time()
    print("Elapsed time:", end - start)

if __name__ == '__main__':
    main(memory_limit)