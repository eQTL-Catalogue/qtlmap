#!/usr/bin/env python

import duckdb
import argparse

def join_files(merged_susie_file, sumstat_batch_file, chrom, start_pos, end_pos, output_file):
    con = duckdb.connect()
    query = f"""
        COPY (
            SELECT a.molecular_trait_id,
              a.chromosome,
              a.position,
              a.ref,
              a.alt,
              a.variant,
              a.ma_samples,
              a.maf,
              a.pvalue,
              a.beta,
              a.se,
              a.type,
              a.ac,
              a.r2,
              a.molecular_trait_object_id,
              a.gene_id,
              a.median_tpm,
              a.rsid,b.cs_id,
              b.cs_size,
              b.pip,
              b.z,
              b.cs_min_r2,
              b.region
            FROM read_parquet('{sumstat_batch_file}') AS a
            INNER JOIN (
                SELECT * FROM read_parquet('{merged_susie_file}')
                WHERE chromosome = '{chrom}' 
                  AND position BETWEEN {start_pos} AND {end_pos}
            ) AS b
            ON a.variant = b.variant AND a.molecular_trait_id = b.molecular_trait_id
        ) TO '{output_file}' (FORMAT 'parquet');
    """
    
    con.execute(query)
    con.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Join filtered run_susie with nominal parquet files using DuckDB.")
    parser.add_argument('-r', '--merged_susie_file', required=True, help="Path to the merged susie parquet file.")
    parser.add_argument('-n', '--sumstat_btach_file', required=True, help="Path to the sumstat batch parquet file.")
    parser.add_argument('-c', '--chrom', required=True, help="Chromosome number of the nominal file.")
    parser.add_argument('-s', '--start_pos', required=True, type=int, help="Start position range of the nominal file.")
    parser.add_argument('-e', '--end_pos', required=True, type=int, help="End position range of the nominal file.")
    parser.add_argument('-o', '--output_file', required=True, help="Path to the output parquet file.")
    
    args = parser.parse_args()
    
    join_files(args.merged_susie_file, args.sumstat_btach_file, args.chrom, args.start_pos, args.end_pos, args.output_file)
