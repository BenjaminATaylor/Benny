#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//where to save the downloaded fastq files and initial versions of bam files (temporary, removed at the end)
params.saveTemp = "/depot/bharpur/data/popgenomes/nextflow"
//where to save the output files generated in the analyses
params.savePath = "/depot/bharpur/data/popgenomes/nextflow/"

//reference genome
params.refGenome = "/depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna"

//NCBI/EBI database taxonomic id for target species
//default is Apis mellifera
params.taxID = "7460"

//path containing the fastq files
//this is correct if you follow the readme guide and download to scratch
//pattern for filenaming
params.fastqPath = "$RCAC_SCRATCH/fastq_files/"
params.fastqPattern = "SRR*_{1,2}.fastq.gz"
/*
========================================================================================
========================================================================================
//A GENERALIZED VERSION OF THE DOWNLOAD, ALIGN, AND CALL HAPLOTYPES PIPELINE
//ORIGINALLY DEVELOPED FOR HARPUR LAB (PURDUE ENTOMOLOGY) FOR APIS MELLIFERA GENOMICS
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
    errorStrategy 'ignore'

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
    errorStrategy 'ignore'

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
    errorStrategy 'ignore'

    input:
    val runAccession

    output:
    val runAccession

    script:
    """
    gatk MarkDuplicates\
        -I ${params.saveTemp}/updated_bam_files/"$runAccession"_updated.bam \
        -O ${params.savePath}/final_bam_files/"$runAccession"_final.bam \
        -M marked_dup_metrics.txt  \
        --REMOVE_SEQUENCING_DUPLICATES true
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
    errorStrategy 'ignore'

    input:
    val runAccession

    output:
    val runAccession

    script:
    """
    gatk BaseRecalibrator \
        -I ${params.savePath}/final_bam_files/"$runAccession"_final.bam \
        -R ${params.refGenome} \
        --known-sites /depot/bharpur/data/ref_genomes/AMELknown_sites.vcf \
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
        --known-sites /depot/bharpur/data/ref_genomes/AMELknown_sites.vcf \
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
    errorStrategy 'ignore'

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
        --known-sites /depot/bharpur/data/ref_genomes/AMELknown_sites.vcf \
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
    errorStrategy 'ignore'

    input:
    val runAccession

    output:
    val runAccession

    script:
    """   
    #THREE
    gatk ApplyBQSR \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam \
        -R ${params.refGenome} \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table \
        -O ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_3.bam
        
    gatk BaseRecalibrator \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_3.bam  \
        -R ${params.refGenome} \
        --known-sites /depot/bharpur/data/ref_genomes/AMELknown_sites.vcf \
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
    errorStrategy 'ignore'

    input:
    val runAccession

    output:
    val runAccession

    script:
    """
    gatk HaplotypeCaller \
        -R ${params.refGenome} \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_3.bam \
        -ERC GVCF \
        -O ${params.savePath}/raw_snps/"$runAccession"_raw_snps.vcf
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

    base_recal1(dupl_removed) \
        | set {base_recal1}
    
    base_recal2(base_recal1) \
        | set {base_recal2}

    base_recal3(base_recal2) \
        | set {base_recal_done}

    haplotype_caller(base_recal_done) \
        | set {haplotype_call}

    alignment_cleanup(haplotype_call) \
        | view()
}