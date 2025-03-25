#!/usr/bin/env python

import argparse
import duckdb


def extract_unique_molecular_trait_ids(concatenated_susie_output, output_file,memory_limit):
    memory_limit = f"{float(memory_limit) * 0.8:.1f}"
    duckdb.sql(f"""
    SET memory_limit='{memory_limit}GB';
    COPY (
        SELECT DISTINCT molecular_trait_id
        FROM read_parquet('{concatenated_susie_output}')
    ) TO '{output_file}' (FORMAT PARQUET)
""")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate pq file with unique molecular trait ids using DuckDB.")
    parser.add_argument('-f', '--concatenated_susie_output', required=True, help="Merged susie output file.")
    parser.add_argument('-o', '--output', required=True, help="Output parquet file.")
    parser.add_argument('-m', '--memory_limit', required=True, type=str, help="Memory limit in GB for DuckDB.")

    args = parser.parse_args()
    extract_unique_molecular_trait_ids(args.concatenated_susie_output, args.output, args.memory_limit)
