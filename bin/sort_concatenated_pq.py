#!/usr/bin/env python

import duckdb
import argparse

def sort_parquet_files(input_file, memory_limit,name_prefix):
    sorted_output_file = name_prefix + '_merged.parquet'
    memory_limit = str(memory_limit*0.8)
    duckdb.sql(f"""
        SET memory_limit='{memory_limit}GB';
        COPY (
            SELECT * FROM read_parquet('{input_file}')
            ORDER BY chromosome, position
        ) TO '{sorted_output_file}' (FORMAT PARQUET)
    """)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple parquet files using DuckDB.")
    parser.add_argument('-i', '--input_file', required=True, help="Input parquet file.")
    parser.add_argument('-m', '--memory_limit', required=True, type=int, help="Memory limit in GB for DuckDB.")
    parser.add_argument('-n', '--name_prefix', required=True, help="Output parquet file prefix.")


    args = parser.parse_args()
    sort_parquet_files(args.input_file,args.memory_limit,args.name_prefix)