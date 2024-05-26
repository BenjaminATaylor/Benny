/*
========================================================================================
    HAPLOTYPE CALLER
========================================================================================
*/

process haplotype_caller{
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 1-00:00:00 -A bharpur'
    //errorstrategy 'ignore'
    
    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(bam), val(runAccession)

    output:
    val runAccession, emit: accession
    val vcfname, emit: gvcf
    
    script:
    vcfname = params.savePath + "/raw_snps/" + runAccession + "_raw_snps.vcf"
    """
    gatk HaplotypeCaller \
        -R $refgenome \
        -I $bam \
        -ERC GVCF \
        -O ${params.savePath}/raw_snps/"$runAccession"_raw_snps.vcf
    """
}

process combine_gvcfs{
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 14 --mem=100G --time 1-00:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    val vcfslist

    output:
    path 'combined.g.vcf'

    script:
    """    
    gatk CombineGVCFs \
        -R $refgenome \
        --variant $vcfslist \
        -O combined.g.vcf
    """
}

process genotype_gvcfs{
    
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 14 --mem=100G --time 1-00:00:00 -A bharpur'
    publishDir "${params.savePath}/raw_snps", mode: 'copy'
    
    input:
    tuple path(refgenome), path(refindex), path(refdict)
    path combinedgvcfs

    output:
    tuple path('combined.vcf.gz'),path('combined.vcf.gz.tbi')
    
    script:
    """
    gatk GenotypeGVCFs \
        -R $refgenome \
        -V $combinedgvcfs \
        -O combined.vcf.gz    
    """

    
}