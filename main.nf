#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//where to save the downloaded fastq files and initial versions of bam files (temporary, removed at the end)
params.saveTemp = "/depot/bharpur/data/popgenomes/nextflow"
//where to save the output files generated in the analyses
params.savePath = "/depot/bharpur/data/popgenomes/nextflow/"

//reference genome (bee RefSeq by default)
params.refGenome = "/depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna"

//set knownsites empty by default
params.knownSites = null

//path containing the fastq files
//this is correct if you follow the readme guide and download to scratch
//pattern for filenaming
params.fastqPath = "$RCAC_SCRATCH/fastq_files/"
params.fastqPattern = "*_{1,2}.fq"
/*
========================================================================================
========================================================================================
//A GENERALIZED VERSION OF THE DOWNLOAD, ALIGN, AND CALL HAPLOTYPES PIPELINE
//ORIGINALLY DEVELOPED FOR HARPUR LAB (PURDUE ENTOMOLOGY) FOR APIS MELLIFERA GENOMICS BY BENNY GOLLER, OCT 2022
//UPDATED BY BENJAMIN A TAYLOR, MAY 2024
========================================================================================
========================================================================================
*/

/*
========================================================================================
========================================================================================
//PROCESSES FOR THE PIPELINE
//PIPELINE EXECUTION CONTROL IN WORKFLOW AT THE END
========================================================================================
    SETUP AND RUN ALIGNMENT
========================================================================================
*/


process alignment{
    tag "$runAccession"
    module 'bioinfo:samtools'
    clusterOptions '--ntasks 1 --time 1-0:00:00 --mem=4G -A bharpur'
    //errorstrategy 'ignore'

    input:
    tuple val(runAccession), file(fastqs)

    output:
    val runAccession

    script:
    """
    mkdir -p ${params.saveTemp}/bam_files
    mkdir -p ${params.saveTemp}/updated_bam_files
    mkdir -p ${params.saveTemp}/final_bam_files
    mkdir -p ${params.saveTemp}/data_tables
    mkdir -p ${params.saveTemp}/recal_bam_files
    mkdir -p ${params.saveTemp}/recal_plots
    mkdir -p ${params.saveTemp}/raw_snps
    /depot/bharpur/apps/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm -r ${params.refGenome} -1 $fastqs[0] -2 $fastqs[1] | samtools sort > ${params.saveTemp}/bam_files/"$runAccession".bam
    """
}
/*
========================================================================================
    DUPLICATE CONTROL
========================================================================================
*/

process check_duplicates{
    tag "$runAccession"
    module 'bioinfo:samtools:picard-tools:GATK'
    clusterOptions '--ntasks 1 --time 4:00:00 -A bharpur'
    //errorstrategy 'ignore'

    input:
    val runAccession

    output:
    val runAccession

    script:
    """
    java -jar \$CLASSPATH AddOrReplaceReadGroups \
        -I ${params.saveTemp}/bam_files/"$runAccession".bam \
        -O ${params.saveTemp}/updated_bam_files/"$runAccession"_updated.bam \
        -RGID 1 \
        -RGLB lib1 \
        -RGPL illumina \
        -RGPU unit1 \
        -RGSM $runAccession           
    """
}

process remove_duplicates{
    tag "$runAccession"
    module 'bioinfo:samtools:picard-tools:GATK'
    clusterOptions '--ntasks 1 --time 4:00:00 -A bharpur'
    //errorstrategy 'ignore'

    input:
    val runAccession

    output:
    tuple path('duprem.bam'),val(runAccession)

    script:
    """
    mkdir -p ${params.saveTemp}/fastq_files
    mkdir -p ${params.saveTemp}/bam_files
    mkdir -p ${params.saveTemp}/updated_bam_files
    mkdir -p ${params.savePath}/final_bam_files
    mkdir -p ${params.savePath}/data_tables
    mkdir -p ${params.saveTemp}/recal_bam_files
    mkdir -p ${params.savePath}/recal_plots
    mkdir -p ${params.savePath}/raw_snps
    
    gatk MarkDuplicates\
        -I ${params.saveTemp}/updated_bam_files/"$runAccession"_updated.bam \
        -O duprem.bam \
        -M marked_dup_metrics.txt  \
        --REMOVE_SEQUENCING_DUPLICATES true \
        --CREATE_INDEX true
    """
}

/*
========================================================================================
    BASE RECALIBRATION
========================================================================================
*/
process base_recal1{
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'
    //errorstrategy 'ignore'

    input:
    tuple path(bam), val(runAccession)

    output:
    val runAccession

    script:
    """        
    gatk BaseRecalibrator \
        -I $bam \
        -R ${params.refGenome} \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table


    #ONE
    gatk ApplyBQSR \
        -I ${params.savePath}/final_bam_files/"$runAccession"_final.bam \
        -R ${params.refGenome} \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table \
        -O ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1.bam

    gatk BaseRecalibrator \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1.bam  \
        -R ${params.refGenome} \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_2.table
    
    gatk AnalyzeCovariates \
        -before ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table \
        -after ${params.savePath}/data_tables/"$runAccession"_recal_data_2.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_1.pdf
    """
}

process base_recal2{
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'
    //errorstrategy 'ignore'

    input:
    val runAccession

    output:
    val runAccession

    script:
    """
    #TWO
    gatk ApplyBQSR \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1.bam \
        -R ${params.refGenome} \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_2.table \
        -O ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam
        
    gatk BaseRecalibrator \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam  \
        -R ${params.refGenome} \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table
    
    gatk AnalyzeCovariates \
        -before ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table \
        -after ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_2.pdf

    rm ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1*
    """
}

process base_recal3{
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'
    //errorstrategy 'ignore'

    input:
    val runAccession

    output:
    tuple path('finalrecal.bam'), val(runAccession)

    script:
    """   
    #THREE
    gatk ApplyBQSR \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam \
        -R ${params.refGenome} \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table \
        -O finalrecal.bam
        
    gatk BaseRecalibrator \
        -I finalrecal.bam  \
        -R ${params.refGenome} \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_4.table
    
    gatk AnalyzeCovariates \
        -before ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table \
        -after ${params.savePath}/data_tables/"$runAccession"_recal_data_4.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_3.pdf

    rm ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2*
    """
}
/*
========================================================================================
    HAPLOTYPE CALLER
========================================================================================
*/
process haplotype_caller{
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 4-00:00:00 -A bharpur'
    //errorstrategy 'ignore'
    
    input:
    tuple path(bam), val(runAccession)

    output:
    val runAccession, emit: accession
    val vcfname, emit: gvcf
    
    script:
    vcfname = params.savePath + "/raw_snps/" + runAccession + "_raw_snps.vcf"
    """
    gatk HaplotypeCaller \
        -R ${params.refGenome} \
        -I $bam \
        -ERC GVCF \
        -O ${params.savePath}/raw_snps/"$runAccession"_raw_snps.vcf
    """
}

process combine_gvcfs{
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 14 --mem=100G --time 4-00:00:00 -A bharpur'

    input:
    val vcfslist

    output:
    path 'combined.g.vcf'

    script:
    """    
    gatk CombineGVCFs \
        -R ${params.refGenome} \
        --variant $vcfslist \
        -O combined.g.vcf
    """
}

process genotype_gvcfs{
    
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 14 --mem=100G --time 4-00:00:00 -A bharpur'
    publishDir "${params.savePath}/raw_snps", mode: 'copy'
    
    input:
    path combinedgvcfs

    output:
    path 'combined.vcf.gz'
    
    script:
    """
    gatk GenotypeGVCFs \
        -R ${params.refGenome} \
        -V $combinedgvcfs \
        -O combined.vcf.gz    
    """

    
}

/*
========================================================================================
    CLEAN-UP
========================================================================================
*/

process alignment_cleanup{
    tag "$runAccession"
    clusterOptions '--ntasks 1 --time 00:10:00 -A bharpur'

    input:
    val runAccession

    output:
    val runAccession
    
    script:
    """
    rm ${params.saveTemp}/bam_files/"$runAccession".bam                        
    rm ${params.saveTemp}/updated_bam_files/"$runAccession"_updated.bam
    rm ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_3*
    """
}
/*
========================================================================================
========================================================================================
PIPELINE EXECUTION CONTROL
*/
workflow{
    Channel.fromFilePairs(params.fastqPath+params.fastqPattern)
        | view() \
        | set {fastq_secured}

    alignment(fastq_secured) \
        | set{align_done}
    
    check_duplicates(align_done) \
        | set {dupl_check_complete}

    remove_duplicates(dupl_check_complete) \
        | set {dupl_removed}

    if( params.knownSites != null ) {
        base_recal1(dupl_removed) \
            | set {base_recal1}
        base_recal2(base_recal1) \
            | set {base_recal2}
        base_recal3(base_recal2) \
            | set {base_recal_done}
        haplotype_caller(base_recal_done) 
    } else {
        haplotype_caller(dupl_removed)        
    }
    
    haplotype_caller.out.gvcf.collectFile(name: 'vcfs.list', newLine: true) \
        | set {vcfslist}
        
    combine_gvcfs(vcfslist)
    genotype_gvcfs(combine_gvcfs.out) 
    
    genotype_gvcfs.out.view()

    //alignment_cleanup(haplotype_caller.out.accession)
}