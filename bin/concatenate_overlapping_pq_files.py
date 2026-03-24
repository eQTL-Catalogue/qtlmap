#!/usr/bin/env python

import duckdb
import argparse

def concatenate_parquet_files(input_files, output_file, memory_limit, threads):
    input_files_list = ', '.join([f"'{file}'" for file in input_files])
    memory_limit = int(memory_limit * 0.9)
    con = duckdb.connect()
    columns_info = con.execute(f"DESCRIBE SELECT * FROM read_parquet('{input_files[0]}')").fetchall()
    column_names = [col[0] for col in columns_info]
    group_keys = ['variant', 'molecular_trait_id', 'rsid']
    select_columns = []
    for col in column_names:
        if col in group_keys:
            select_columns.append(col)
        else:
            select_columns.append(f"ANY_VALUE({col}) AS {col}")
    column_syntax_to_select = ", ".join(select_columns)
    print(column_syntax_to_select)
    con.execute("PRAGMA enable_profiling='json';")
    con.execute(f"SET memory_limit='{int(memory_limit)}MB'")
    con.execute(f"PRAGMA threads={threads}")
    con.execute("PRAGMA profiling_output='profile.json';")
    con.execute(f"""
        COPY (
         SELECT {column_syntax_to_select}
         FROM read_parquet([{input_files_list}])
        GROUP BY variant, molecular_trait_id, rsid
            ORDER BY chromosome, position
        ) TO '{output_file}' (FORMAT PARQUET);
    """)

    con.close()
    print(f"Merged and sorted files into {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Concatenate multiple parquet files using DuckDB.")
    parser.add_argument('-f', '--files', nargs='+', required=True, help="List of parquet files to concatenate.")
    parser.add_argument('-o', '--output', required=True, help="Output parquet file.")
    parser.add_argument('-m', '--memory_limit', required=True, type=int, help="Memory limit in MB for DuckDB.")
    parser.add_argument('-t', '--threads', type=int, default=2, help="Number of CPU threads for DuckDB (default: 2)")
    args = parser.parse_args()
    concatenate_parquet_files(args.files, args.output, args.memory_limit, args.threads)