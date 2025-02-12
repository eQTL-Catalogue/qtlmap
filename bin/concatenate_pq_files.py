#!/usr/bin/env python

import duckdb
import argparse

def concatenate_parquet_files(input_files, output_file,memory_limit):
    input_files_list = ', '.join([f"'{file}'" for file in input_files])
    memory_limit = f"{float(memory_limit) * 0.8:.1f}"
    query = f"""
        SET memory_limit='{memory_limit}GB';
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
    parser.add_argument('-m', '--memory_limit', required=True, type=str, help="Memory limit in GB for DuckDB.")


    args = parser.parse_args()
    concatenate_parquet_files(args.files, args.output,args.memory_limit)