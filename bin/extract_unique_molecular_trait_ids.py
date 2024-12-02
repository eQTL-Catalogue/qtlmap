#!/usr/bin/env python

import argparse
import duckdb


def extract_unique_molecular_trait_ids(concatenated_susie_output, output_file):
    duckdb.sql(f"""
    COPY (
        SELECT DISTINCT molecular_trait_id
        FROM read_parquet('{concatenated_susie_output}')
    ) TO '{output_file}' (FORMAT PARQUET)
""")
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate pq file with unique molecular trait ids using DuckDB.")
    parser.add_argument('-f', '--concatenated_susie_output', required=True, help="Merged susie output file.")
    parser.add_argument('-o', '--output', required=True, help="Output parquet file.")

    args = parser.parse_args()
    extract_unique_molecular_trait_ids(args.concatenated_susie_output, args.output)
