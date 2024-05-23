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

include { haplotype_caller; haplotype_caller as haplotype_caller_2 } from './modules.nf'
include { combine_gvcfs; combine_gvcfs as combine_gvcfs_2 } from './modules.nf'
include { genotype_gvcfs; genotype_gvcfs as genotype_gvcfs_2 } from './modules.nf'

ch_refgenome = Channel.value(file(params.refGenome))

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
    clusterOptions '--ntasks 1 --time 08:00:00 --mem=4G -A bharpur'
    //errorstrategy 'ignore'

    input:
    path refgenome
    tuple val(runAccession), file(fastqs)

    output:
    tuple val(runAccession), path('init.bam')

    script:
    """
    mkdir -p ${params.saveTemp}/bam_files
    mkdir -p ${params.saveTemp}/updated_bam_files
    mkdir -p ${params.saveTemp}/final_bam_files
    mkdir -p ${params.saveTemp}/data_tables
    mkdir -p ${params.saveTemp}/recal_bam_files
    mkdir -p ${params.saveTemp}/recal_plots
    mkdir -p ${params.saveTemp}/raw_snps
    /depot/bharpur/apps/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm -r $refgenome -1 $fastqs[0] -2 $fastqs[1] | samtools sort > init.bam
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
    tuple val(runAccession), path(inbam)

    output:
    tuple val(runAccession), path('AORRG.bam')

    script:
    """
    java -jar \$CLASSPATH AddOrReplaceReadGroups \
        -I $inbam \
        -O AORRG.bam \
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
    tuple val(runAccession), path(inbam)

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
        -I $inbam \
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
    path refgenome
    tuple path(bam), val(runAccession)

    output:
    val runAccession

    script:
    """        
    gatk BaseRecalibrator \
        -I $bam \
        -R $refgenome \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table


    #ONE
    gatk ApplyBQSR \
        -I ${params.savePath}/final_bam_files/"$runAccession"_final.bam \
        -R $refgenome \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_1.table \
        -O ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1.bam

    gatk BaseRecalibrator \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1.bam  \
        -R $refgenome \
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
    path refgenome
    val runAccession

    output:
    val runAccession

    script:
    """
    #TWO
    gatk ApplyBQSR \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_1.bam \
        -R $refgenome \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_2.table \
        -O ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam
        
    gatk BaseRecalibrator \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam  \
        -R $refgenome \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table
    
    gatk AnalyzeCovariates \
        -before ${params.savePath}/data_tables/"$runAccession"_recal_data_2.table \
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
    path refgenome
    val runAccession

    output:
    tuple path('finalrecal.bam'), val(runAccession)

    script:
    """   
    #THREE
    gatk ApplyBQSR \
        -I ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2.bam \
        -R $refgenome \
        --bqsr-recal-file ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table \
        -O finalrecal.bam
        
    gatk BaseRecalibrator \
        -I finalrecal.bam  \
        -R $refgenome \
        --known-sites ${params.knownSites} \
        -O ${params.savePath}/data_tables/"$runAccession"_recal_data_4.table
    
    gatk AnalyzeCovariates \
        -before ${params.savePath}/data_tables/"$runAccession"_recal_data_3.table \
        -after ${params.savePath}/data_tables/"$runAccession"_recal_data_4.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_3.pdf

    rm ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_2*
    """
}

/*
========================================================================================
    BQSR BOOTSTRAPPING
    (parts of this taken from gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)
========================================================================================
*/

process select_snps{
    
    module 'bioinfo:GATK'
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    path refgenome
    tuple path(rawvcf),path(idx)
    
    output:
    path 'raw_snps.vcf'
    
    script:
    """
    gatk SelectVariants \
        -R $refgenome \
        -V $rawvcf \
        -select-type SNP \
        -O raw_snps.vcf
    """ 
}

process select_indels{
    
    module 'bioinfo:GATK'
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    path refgenome
    tuple path(rawvcf),path(idx)
    
    output:
    path 'raw_indels.vcf'
    
    script:
    """
    gatk SelectVariants \
        -R $refgenome \
        -V $rawvcf \
        -select-type INDEL \
        -O raw_indels.vcf
    """ 
}

process filter_snps{
    
    module 'bioinfo:GATK'
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    path refgenome
    path rawsnpsvcf
    
    output:
    path 'filtered_snps.vcf'
    
    script:
    """
    #Note that here we are looking only for the highest-confidence SNPs for downstream filtering, so we're going to use more conservative filters than typical
    gatk VariantFiltration \
        -R $refgenome \
        -V $rawsnpsvcf \
        -O filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 25.0" \
    """ 
}

process filter_indels{
    
    module 'bioinfo:GATK'
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    path refgenome
    path rawindelsvcf
    
    output:
    path 'filtered_indels.vcf'
    
    script:
    """
    #Note that here we are looking only for the highest-confidence InDels for downstream filtering, so we're going to use more conservative filters than typical
    gatk VariantFiltration \
        -R $refgenome \
        -V $rawindelsvcf \
        -O filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 25.0" \
    """ 
}

process apply_snp_filter{
    
    module 'bioinfo:GATK'
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    path filteredsnps
    
    output:
    tuple path('snps_filtered.vcf'), path('snps_filtered.vcf.idx')
    
    script:
    """
    gatk SelectVariants \
        --exclude-filtered \
        -V $filteredsnps \
        -O snps_filtered.vcf
    """
}

process apply_indel_filter{
    
    module 'bioinfo:GATK'
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    path rawindelsvcf
    
    output:
    tuple path('indels_filtered.vcf'), path('indels_filtered.vcf.idx')
    
    script:
    """
    gatk SelectVariants \
        --exclude-filtered \
        -V $rawindelsvcf \
        -O indels_filtered.vcf
    """
}

process base_recal_boot_init{
    
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'

    input:
    path refgenome 
    tuple path(bam), val(runAccession)
    tuple path(snps), path(snp_idx)
    tuple path(indels), path(indel_idx)

    output:
    tuple path(bam), val(runAccession), path('recal_data.table'), path(snps), path(snp_idx), path(indels), path(indel_idx)
    
    script:
    """
	gatk BaseRecalibrator \
        -R $refgenome \
        -I $bam \
        --known-sites $snps \
        --known-sites $indels \
        -O recal_data.table
    """
} 

       
process base_recal_boot_1{
    
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'

    input:
    path refgenome
    tuple path(bam), val(runAccession), path(table), path(snps), path(snp_idx), path(indels), path(indel_idx)


    output:
    tuple path('recal.bam'), val(runAccession), path('recal_data.table'), path(snps), path(snp_idx), path(indels), path(indel_idx)
    
    script:        
    """
    #ONE
    gatk ApplyBQSR \
        -I $bam \
        -R $refgenome \
        --bqsr-recal-file $table \
        -O recal.bam

    gatk BaseRecalibrator \
        -I recal.bam  \
        -R $refgenome \
        --known-sites $snps \
        --known-sites $indels \
        -O recal_data.table
            
    gatk AnalyzeCovariates \
        -before $table \
        -after recal_data.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_1.pdf
    """
} 

process base_recal_boot_2{
    
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'

    input:
    path refgenome
    tuple path(bam), val(runAccession), path(table), path(snps), path(snp_idx), path(indels), path(indel_idx)

    output:
    tuple path('recal_2.bam'), val(runAccession), path('recal_data_2.table'), path(snps), path(snp_idx), path(indels), path(indel_idx)
    
    script:        
    """
    #TWO
    gatk ApplyBQSR \
        -I $bam \
        -R $refgenome \
        --bqsr-recal-file $table \
        -O recal_2.bam

    gatk BaseRecalibrator \
        -I recal_2.bam  \
        -R $refgenome \
        --known-sites $snps \
        --known-sites $indels \
        -O recal_data_2.table
            
    gatk AnalyzeCovariates \
        -before $table \
        -after recal_data_2.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_2.pdf
    """
} 

process base_recal_boot_3{
    
    tag "$runAccession"
    module 'bioinfo:GATK'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'

    input:
    path refgenome
    tuple path(bam), val(runAccession), path(table), path(snps), path(snp_idx), path(indels), path(indel_idx)


    output:
    tuple path('recal_3.bam'), val(runAccession), path('recal_data_3.table'), path(snps), path(snp_idx), path(indels), path(indel_idx)
    
    script:        
    """
    #THREE
    gatk ApplyBQSR \
        -I $bam \
        -R $refgenome \
        --bqsr-recal-file $table \
        -O recal_3.bam

    gatk BaseRecalibrator \
        -I recal_3.bam  \
        -R $refgenome \
        --known-sites $snps \
        --known-sites $indels \
        -O recal_data_3.table
            
    gatk AnalyzeCovariates \
        -before $table \
        -after recal_data_3.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_3.pdf
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

    alignment(ch_refgenome,fastq_secured) \
        | set{align_done}
    
    check_duplicates(align_done) \
        | set {dupl_check_complete}

    remove_duplicates(dupl_check_complete) \
        | set {dupl_removed}

    if( params.knownSites != null ) {
        base_recal1(ch_refgenome, dupl_removed) \
            | set {base_recal1}
        base_recal2(ch_refgenome, base_recal1) \
            | set {base_recal2}
        base_recal3(ch_refgenome, base_recal2) \
            | set {base_recal_done}
        haplotype_caller(ch_refgenome, base_recal_done) 
    } else {
        haplotype_caller(ch_refgenome, dupl_removed)        
    }
    
    haplotype_caller.out.gvcf.collectFile(name: 'vcfs.list', newLine: true) \
        | set {vcfslist}
        
    combine_gvcfs(ch_refgenome, vcfslist)
    genotype_gvcfs(ch_refgenome, combine_gvcfs.out) 
    
    // If no known-sites file is provided, we apply a harsh variant filter and use those high-confidence sites to perform base quality recalibration
    if( params.knownSites == null ) {
        
        select_snps(ch_refgenome, genotype_gvcfs.out)
        filter_snps(ch_refgenome, select_snps.out)
        apply_snp_filter(filter_snps.out)

        select_indels(ch_refgenome, genotype_gvcfs.out)
        filter_indels(ch_refgenome, select_indels.out)
        apply_indel_filter(filter_indels.out)
            
        // Note that we use the first() operator here to ensure that the snp/indel filter outputs are treated as values instead of queues
        base_recal_boot_init(
            ch_refgenome, 
            dupl_removed,
            apply_snp_filter.out.first(),
            apply_indel_filter.out.first())
            
        base_recal_boot_1(ch_refgenome, base_recal_boot_init.out)
        base_recal_boot_2(ch_refgenome, base_recal_boot_1.out)
        base_recal_boot_3(ch_refgenome, base_recal_boot_2.out)
        base_recal_boot_3.out.view()
        
        haplotype_caller_2(ch_refgenome, base_recal_boot_3.out) 
        haplotype_caller_2.out.gvcf.collectFile(name: 'vcfs.list', newLine: true) \
        | set {vcfslist_2}
        
        combine_gvcfs_2(ch_refgenome, vcfslist_2)
        genotype_gvcfs_2(ch_refgenome, combine_gvcfs_2.out)
        
        genotype_gvcfs_2.out.view()
    }
    
    //alignment_cleanup(haplotype_caller.out.accession)
}