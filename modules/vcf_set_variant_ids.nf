process vcf_set_variant_ids {
    tag "${qtl_subset}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file(vcf)

    output:
    tuple val(qtl_subset), file("${vcf.simpleName}_renamed.vcf.gz")

    script:
    """
    # originally:  bcftools annotate --set-id 'chr%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' $vcf -Oz -o ${vcf.simpleName}_renamed.vcf.gz
    
    # to make sure that the chromosomes are normalized
    cat > chr_map.txt <<'EOF'
    chr1 1
    chr2 2
    chr3 3
    chr4 4
    chr5 5
    chr6 6
    chr7 7
    chr8 8
    chr9 9
    chr10 10
    chr11 11
    chr12 12
    chr13 13
    chr14 14
    chr15 15
    chr16 16
    chr17 17
    chr18 18
    chr19 19
    chr20 20
    chr21 21
    chr22 22
    chrX X
    chrY Y
    chrM M
    EOF
    
    bcftools annotate --rename-chrs chr_map.txt $vcf -Ou | bcftools annotate --set-id 'chr%CHROM\\_%POS\\_%REF\\_%FIRST_ALT' -Oz -o ${vcf.simpleName}_renamed.vcf.gz
    """
}