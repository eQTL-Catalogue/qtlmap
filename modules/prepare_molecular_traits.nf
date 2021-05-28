/*
 * Prepare molecular trait data by:
 * 1. Filter molecular traits with less than mincisvaraint variants in cis
 * 2. Perform PCA on the molecular trait matrix
 */
process prepare_molecular_traits {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(expression_matrix), file(phenotype_metadata), file(sample_metadata), file(vcf_variant_info)

    output: 
    tuple val(qtl_subset), file("*.bed"), emit: bed_file
    tuple val(qtl_subset), file("*.sample_names.txt"), emit: sample_names
    tuple val(qtl_subset), file("${qtl_subset}.pheno_cov.txt"), emit: pheno_cov
    tuple val(qtl_subset), file(phenotype_metadata), emit: pheno_meta

    script:
    """
    Rscript $baseDir/bin/prepare_molecular_traits.R \\
        -p "$phenotype_metadata" \\
        -s "$sample_metadata" \\
        -e "$expression_matrix" \\
        -v "$vcf_variant_info" \\
        -o "." \\
        -c ${params.cis_window} \\
        -m ${params.mincisvariant} \\
        -a ${params.covariates}
    
    #Merge phenotype covariates together
    head -n ${params.n_pheno_pcs + 1} phenoPCA.tsv > ${qtl_subset}.pheno_cov.txt
    tail -n+2 additional_covariates.tsv >> ${qtl_subset}.pheno_cov.txt
    """
}

// Compress and index to original bed file and make one that is comaptible with fastQTL
process compress_bed {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(bed_file)

    output:
    tuple val(qtl_subset), file("${bed_file}.gz"), file("${bed_file}.gz.tbi"), file("${bed_file.baseName}.fastQTL.bed.gz"), file("${bed_file.baseName}.fastQTL.bed.gz.tbi")

    script:
    """
    bgzip ${bed_file} && tabix -p bed ${bed_file}.gz
    
    csvtk cut -C\$ -t -f -strand,-group_id ${bed_file}.gz | bgzip > ${bed_file.baseName}.fastQTL.bed.gz
    tabix -p bed ${bed_file.baseName}.fastQTL.bed.gz
    """
}

/*
 * STEP 5 - Perform PCA on the genotype and phenotype data
 */
process make_pca_covariates {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/PCA/${qtl_subset}", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(phenotype_cov), file(vcf)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.covariates.txt")

    script:
    """
    plink2 --vcf $vcf --vcf-half-call h --indep-pairwise 50000 200 0.05 --out ${qtl_subset}_pruned_variants --threads ${task.cpus} --memory ${task.memory.mega} --const-fid 
    plink2 --vcf $vcf --vcf-half-call h --extract ${qtl_subset}_pruned_variants.prune.in --make-bed --out ${qtl_subset}_pruned --const-fid 
    plink2 -bfile ${qtl_subset}_pruned --pca ${params.n_geno_pcs} header tabs
    cat plink.eigenvec \\
        | sed '1s/IID/genotype_id/' \\
        | sed '1s/PC/geno_PC/g' \\
        | csvtk cut -t -f -"FID" \\
        | csvtk transpose -t > ${qtl_subset}.geno.pca
    cat $phenotype_cov > ${qtl_subset}.covariates.txt    
    set +o pipefail; tail -n+2 ${qtl_subset}.geno.pca | head -n ${params.n_geno_pcs} >> ${qtl_subset}.covariates.txt
    """
}