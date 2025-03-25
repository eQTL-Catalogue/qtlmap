#!/usr/bin/env python

import argparse
import duckdb


def extract_lead_cc_signal(qtl_ss, unique_ids_pq, output_file_name,memory_limit):
    memory_limit = f"{float(memory_limit) * 0.8:.1f}"
    duckdb.sql(f"""
    SET memory_limit='{memory_limit}GB';
    COPY (
        SELECT qtl_ss.*
        FROM read_parquet('{qtl_ss}') qtl_ss
        INNER JOIN read_parquet('{unique_ids_pq}') unique_ids
        ON qtl_ss.molecular_trait_id = unique_ids.molecular_trait_id
    ) TO '{output_file_name}' (FORMAT PARQUET)
""")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate cc files.")
    parser.add_argument('-s', '--qtl_nom_sum', required=True, help="Nominal summary statistics file")
    parser.add_argument('-u', '--unique_molecular_trait_ids', required=True, help="Parquet file with lead signal (phenotype id) generated by build_connected_components process")
    parser.add_argument('-o', '--output', required=True, help="Output parquet file.")
    parser.add_argument('-m', '--memory_limit', required=True, type=str, help="Memory limit in GB for DuckDB.")


    args = parser.parse_args()
    extract_lead_cc_signal(args.qtl_nom_sum, args.unique_molecular_trait_ids, args.output,args.memory_limit) 
