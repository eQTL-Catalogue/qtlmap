#!/usr/bin/env python
import time
import gzip
import csv
import io
import os
import argparse
import pandas as pd

def make_var_info_rsid_dict(csv_file_path):
    df = pd.read_csv(csv_file_path, sep='\t', comment='#', header=None, dtype=str, names=['variant', 'type','r2'], usecols=[2, 5, 9], verbose=True, keep_default_na = False)
    df = df.assign(compact_info = df["type"].astype(str) + "," + df["r2"].astype(str)) 
    res_dict = df.set_index('variant').to_dict()['compact_info']
    return res_dict

def make_var_info_dict(csv_file_path):
    df = pd.read_csv(csv_file_path, sep='\t', comment='#', header=None, dtype=str, names=['variant', 'type','ac','an','r2'], usecols=[2, 5, 6, 7, 9], verbose=True, keep_default_na = False)
    df = df.assign(compact_info = df["type"].astype(str) + "#" + df["ac"].astype(str) + "#" + df["an"].astype(str) + "#" + df["r2"].astype(str)) 
    res_dict = df.set_index('variant').to_dict()['compact_info']
    # res_dict["variant"] = "type,ac,an,r2"
    return res_dict

def make_var_rsid_dict(csv_file_path):
    df = pd.read_csv(csv_file_path, sep='\t', header=None, dtype=str, names=['variant', 'rsid'], usecols=[0, 1], verbose=True, keep_default_na = False)
    df['rsid'] = df.groupby(['variant'])['rsid'].transform(lambda x: '#'.join(x))
    df.drop_duplicates()
    df = df.assign(compact_info = df["rsid"].astype(str))
    res_dict = df.set_index('variant').to_dict()['compact_info']
    # res_dict["variant"] = "rsid"
    return res_dict

def make_pheno_meta_dict(csv_file_path):
    df = pd.read_csv(csv_file_path, sep='\t', dtype=str, usecols=["phenotype_id","group_id","gene_id"], verbose=True)
    df = df.assign(compact_info = df["group_id"].astype(str) + "#" + df["gene_id"].astype(str)) 
    res_dict = df.set_index('phenotype_id').to_dict()['compact_info']
    # res_dict["molecular_trait_id"] = "molecular_trait_object_id,gene_id"
    return res_dict

def make_median_tpm_dict(csv_file_path):
    df = pd.read_csv(csv_file_path, sep='\t', dtype=str, usecols=["phenotype_id", "median_tpm"], verbose=True)
    res_dict = df.set_index('phenotype_id').to_dict()['median_tpm']
    # res_dict["molecular_trait_id"] = "median_tpm"
    return res_dict

def main():
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-v', help='The variant info file', required=True)
    argparser.add_argument('-r', help='The variant_id rsid map file', required=False)
    argparser.add_argument('-s', help='The summary_statistics file', required=False)
    argparser.add_argument('-p', help='The phenotype metadata file', required=False)
    argparser.add_argument('-m', help='The median_tpm file', required=False)
    argparser.add_argument('-o', help='The output file name', required=True)
    
    args = argparser.parse_args()
    rsid_map = args.r
    var_info = args.v
    summ_stats = args.s
    pheno_meta = args.p
    median_tpm = args.m
    outfile = args.o  

    i=1
    tic = time.process_time()  
    f_w = gzip.GzipFile(outfile, "wb")
    writer = csv.writer(io.TextIOWrapper(f_w, newline="", write_through=True), delimiter='\t')
    if pheno_meta and summ_stats:
        print("Build the variant information dictionary...")
        tic = time.process_time()  
        var_info_dict = make_var_info_dict(var_info)
        print(dict(list(var_info_dict.items())[0:5]))  
        toc = time.process_time()  
        print("time of building var_info dict: ", toc - tic)

        print("Building the var_rsid dictionary...")
        tic = time.process_time()  
        var_rsid_dict = make_var_rsid_dict(rsid_map)
        print(dict(list(var_rsid_dict.items())[0:1]))  
        toc = time.process_time()  
        print("time of building var_rsid dict: ", toc - tic)

        print("Building phenotype metadata dictionary...")
        tic = time.process_time()  
        pheno_meta_dict = make_pheno_meta_dict(pheno_meta)
        print(dict(list(pheno_meta_dict.items())[0:3]))  
        toc = time.process_time()  
        print("time of building var_rsid dict: ", toc - tic)

        if median_tpm!="null.txt":
            print("Building median_tpm dictionary...")
            tic = time.process_time()  
            median_tpm_dict = make_median_tpm_dict(median_tpm)
            print(dict(list(median_tpm_dict.items())[0:3]))  
            toc = time.process_time()  
            print("time of building var_rsid dict: ", toc - tic)

        with gzip.open(summ_stats, "rt") as f:
            col_names = ["molecular_trait_id","chromosome","position","ref","alt","variant","ma_samples","maf","pvalue","beta","se","type","ac","an","r2","molecular_trait_object_id","gene_id","median_tpm","rsid"]
            writer.writerow(col_names)
            for line in f:
                if i%10000000 == 0:
                    print(i)
                i+=1
                k = line.rstrip().split("\t")
                line_wr = []
                # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,ma_count,maf,pvalue,beta,se

                #Remove ma_count column 
                k.pop(7)
                # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se

                if k[5] in var_info_dict:
                    k.extend(var_info_dict[k[5]].split("#")) #join ac,an,type, and r2 values
                # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se,ac,an,type,r2
                
                if k[0] in pheno_meta_dict:
                    k.extend(pheno_meta_dict[k[0]].split("#")) # join pheno_metadata
                else:
                    k.extend(["NA","NA"]) # join NA's for unfound phenotype id
                # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se,ac,an,type,r2,molecular_trait_object_id,gene_id
                
                if median_tpm!="null.txt" and k[16] in median_tpm_dict:
                    k.append(median_tpm_dict[k[16]]) # join median_tpm
                else:
                    k.append("NA") # join NA's for unfound median_tpm phenotype_id
                # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se,ac,an,type,r2,molecular_trait_object_id,gene_id,median_tpm

                if k[5] in var_rsid_dict:
                    k.extend(var_rsid_dict[k[5]].split(",")[::-1]) # join rsid
                    if "#" in k[-1]: # check if there are multiple rsids
                        rsids = k[-1].strip().split("#")
                        for rsid in rsids:
                            split_list = k[:-1]
                            split_list.append(rsid)
                            writer.writerow(split_list)   
                    else:     
                        writer.writerow(k)
                else:
                    k.extend(["NA"]) # add NA's for missing
                    writer.writerow(k)
                # molecular_trait_id,chromosome,position,ref,alt,variant,ma_samples,maf,pvalue,beta,se,ac,an,type,r2,molecular_trait_object_id,gene_id,median_tpm,rsid
    else:
        print("Building the dictionary")
        tic = time.process_time()  
        var_dict = make_var_info_rsid_dict(var_info)
        print(dict(list(var_dict.items())[0:3]))  
        toc = time.process_time()  
        print("time of building dict: ", toc - tic)
        print("joining the files...")
        with gzip.open(rsid_map, "rt") as f:
            for line in f:
                if i%100000000 == 0:
                    print(i)
                i+=1
                k = line.rstrip().split("\t")
                if k[0] in var_dict:
                    writer.writerow(k)
                
    f_w.close()
    toc = time.process_time()  
    print("time of joining matrices: ", toc - tic)
    
if __name__ == '__main__':
    main()        
