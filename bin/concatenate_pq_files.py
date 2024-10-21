#!/usr/bin/env python

import duckdb
import argparse

def concatenate_parquet_files(input_files, output_file):
    input_files_list = ', '.join([f"'{file}'" for file in input_files])
    query = f"""
        COPY (
            SELECT * FROM read_parquet([{input_files_list}])
            ORDER BY chromosome, position
        ) TO '{output_file}' (FORMAT PARQUET);
    """ 
    con = duckdb.connect()
    con.execute(query)
    print(f"Merged files: {input_files} into {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple parquet files using DuckDB.")
    parser.add_argument('-f', '--files', nargs='+', required=True, help="List of parquet files to concatenate.")
    parser.add_argument('-o', '--output', required=True, help="Output parquet file.")

    args = parser.parse_args()
    concatenate_parquet_files(args.files, args.output)