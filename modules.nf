/*
========================================================================================
    HAPLOTYPE CALLER
========================================================================================
*/

process haplotype_caller{
    
    tag "$runAccession"
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--ntasks 1 --time 1-00:00:00 -A bharpur'
    
    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(bam), val(runAccession)

    output:
    val runAccession, emit: accession
    val vcfname, emit: gvcf
    
    script:
    vcfname = params.savePath + "/raw_snps/" + runAccession + "_raw_snps.vcf.gz"
    """
    gatk HaplotypeCaller \
        -R $refgenome \
        -I $bam \
        -ERC GVCF \
        -O ${params.savePath}/raw_snps/"$runAccession"_raw_snps.vcf.gz
    """
}

process combine_gvcfs{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--ntasks 14 --mem=100G --time 1-00:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    val vcfslist

    output:
    tuple path('combined.g.vcf.gz'), path('combined.g.vcf.gz.tbi')


    script:
    """    
    gatk CombineGVCFs \
        -R $refgenome \
        --variant $vcfslist \
        -O combined.g.vcf.gz
        
    gatk IndexFeatureFile \
        -I combined.g.vcf.gz
    """
}

process genotype_gvcfs{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--ntasks 14 --mem=100G --time 1-00:00:00 -A bharpur'
    publishDir "${params.savePath}/raw_snps", mode: 'copy'
    
    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(combinedgvcfs),path(index)

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