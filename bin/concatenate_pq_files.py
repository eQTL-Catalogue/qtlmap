#!/usr/bin/env python

import duckdb
import argparse

def concatenate_parquet_files(input_files, output_file,memory_limit,threads):
    input_files_list = ', '.join([f"'{file}'" for file in input_files])
    memory_limit = int(memory_limit * 0.9)
    query = f"""
        SET memory_limit='{memory_limit}GB';
        COPY (
            SELECT * FROM read_parquet([{input_files_list}])
            ORDER BY chromosome, position
        ) TO '{output_file}' (FORMAT PARQUET);
    """ 
    con = duckdb.connect()
    con.execute("PRAGMA enable_profiling='json';")
    con.execute(f"SET memory_limit='{int(memory_limit)}MB'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA profiling_output='profile.json';")
    con.execute(query)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple parquet files using DuckDB.")
    parser.add_argument('-f', '--files', nargs='+', required=True, help="List of parquet files to concatenate.")
    parser.add_argument('-o', '--output', required=True, help="Output parquet file.")
    parser.add_argument('-m', '--memory_limit', required=True, type=int, help="Memory limit in MB for DuckDB.")
    parser.add_argument('-t', '--threads', type=int, default=2, help="Number of CPU threads for DuckDB (default: 2)")

    args = parser.parse_args()
    concatenate_parquet_files(args.files, args.output,args.memory_limit,args.threads)