#!/usr/bin/env nextflow

//where to save the output files generated in the analyses
params.savePath = "/depot/bharpur/data/popgenomes/nextflow/"

//reference genome (bee RefSeq by default)
params.refGenome = "/depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna"

//set knownsites empty by default
params.knownSites = null

//path containing the fastq files
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
//HEAVILY UPDATED BY BENJAMIN A TAYLOR, MAY 2024
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

process index_genome{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 5:00 --mem=1G -A bharpur'
    
    input:
    path refgenome
    
    output:
    tuple path(refgenome), path('*.fai'), path('*.dict')
    
    script:
    """
    samtools faidx $refgenome > "$refgenome".fai
    foo=\$(sed -r "s/\\.f[n]a.*//g" <<< $refgenome)
    gatk CreateSequenceDictionary -R $refgenome -O "\$foo.dict"
    """
    
}

process alignment{
    
    tag "$runAccession"
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 04:00:00 --mem=4G -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple val(runAccession), file(fastqs)

    output:
    tuple val(runAccession), path('init.bam')

    script:
    """
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
    //module 'biocontainers:samtools:picard:gatk4'
    container "quay.io/biocontainers/picard:3.1.1--hdfd78af_0"
    clusterOptions '--time 4:00:00 -A bharpur'

    input:
    tuple val(runAccession), path(inbam)

    output:
    tuple val(runAccession), path('AORRG.bam')

    script:
    """
    picard AddOrReplaceReadGroups \
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
    //module 'biocontainers:samtools:picard:gatk4'
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 4:00:00 --mem=4G -A bharpur'

    input:
    tuple val(runAccession), path(inbam)

    output:
    tuple path('duprem.bam'),val(runAccession)

    script:
    """
    mkdir -p ${params.savePath}/final_bam_files
    mkdir -p ${params.savePath}/data_tables
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
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    path knownsites
    tuple path(bam), val(runAccession)

    output:
    tuple val(runAccession), path('recal_1.bam'), path('recal_data_2.table')

    script:
    """        
    gatk BaseRecalibrator \
        -I $bam \
        -R $refgenome \
        --known-sites $knownsites \
        -O recal_data_1.table


    #ONE
    gatk ApplyBQSR \
        -I ${params.savePath}/final_bam_files/"$runAccession"_final.bam \
        -R $refgenome \
        --bqsr-recal-file recal_data_1.table \
        -O recal_1.bam

    gatk BaseRecalibrator \
        -I recal_1.bam  \
        -R $refgenome \
        --known-sites $knownsites \
        -O recal_data_2.table
    
    gatk AnalyzeCovariates \
        -before recal_data_1.table \
        -after $recal_data_2.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_1.pdf
    """
}

process base_recal2{
    
    tag "$runAccession"
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    path knownsites
    tuple val(runAccession), path(recal_1_bam), path(recal_table_2)

    output:
    tuple val(runAccession), path('recal_2.bam'), path('recal_data_3.table')

    script:
    """
    #TWO
    gatk ApplyBQSR \
        -I $recal_1_bam \
        -R $refgenome \
        --bqsr-recal-file $recal_table_2 \
        -O recal_2.bam
        
    gatk BaseRecalibrator \
        -I recal_2.bam  \
        -R $refgenome \
        --known-sites $knownsites \
        -O recal_data_3.table
    
    gatk AnalyzeCovariates \
        -before $recal_table_2 \
        -after recal_data_3.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_2.pdf
    """
}

process base_recal3{
    
    tag "$runAccession"
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'


    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple val(runAccession), path(recal_2_bam), path(recal_table_3)

    output:
    tuple path('finalrecal.bam'), val(runAccession)

    script:
    """   
    #THREE
    gatk ApplyBQSR \
        -I $recal_2_bam \
        -R $refgenome \
        --bqsr-recal-file $recal_table_3 \
        -O finalrecal.bam
        
    gatk BaseRecalibrator \
        -I finalrecal.bam  \
        -R $refgenome \
        --known-sites $knownsites \
        -O recal_data_4.table
    
    gatk AnalyzeCovariates \
        -before $recal_table_3 \
        -after recal_data_4.table \
        -plots ${params.savePath}/recal_plots/"$runAccession"_plot_3.pdf
    """
}

/*
========================================================================================
    BQSR BOOTSTRAPPING
    (parts of this taken from gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/)
========================================================================================
*/

process select_snps{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(rawvcf),path(idx)
    
    output:
    tuple path('raw_snps.vcf.gz'), path('raw_snps.vcf.gz.tbi')
    
    script:
    """
    gatk SelectVariants \
        -R $refgenome \
        -V $rawvcf \
        -select-type SNP \
        -O raw_snps.vcf.gz
        
    gatk IndexFeatureFile \
        -I raw_snps.vcf.gz
    """ 
}

process select_indels{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(rawvcf),path(idx)
    
    output:
    tuple path('raw_indels.vcf.gz'), path('raw_indels.vcf.gz.tbi')

    script:
    """
    gatk SelectVariants \
        -R $refgenome \
        -V $rawvcf \
        -select-type INDEL \
        -O raw_indels.vcf.gz
        
    gatk IndexFeatureFile \
        -I raw_indels.vcf.gz
    """ 
}

process filter_snps{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(rawsnpsvcf),path(index)
    
    output:
    tuple path('snps_filtered.vcf.gz'), path('snps_filtered.vcf.gz.tbi')
    
    script:
    """
    #Note that here we are looking only for the highest-confidence SNPs for downstream filtering, so we're going to use more conservative filters than typical
    gatk VariantFiltration \
        -R $refgenome \
        -V $rawsnpsvcf \
        -O filtered_snps.vcf.gz \
        -filter-name "QD_filter" -filter "QD < 25.0" 
        
    gatk SelectVariants \
        --exclude-filtered \
        -V filtered_snps.vcf.gz \
        -O snps_filtered.vcf.gz
        
    gatk IndexFeatureFile \
        -I snps_filtered.vcf.gz
    """ 
}

process filter_indels{
    
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--mem=50G --time 1-00:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(rawindelsvcf),path(index)
    
    output:
    tuple path('indels_filtered.vcf.gz'), path('indels_filtered.vcf.gz.tbi')
    
    script:
    """
    #Note that here we are looking only for the highest-confidence InDels for downstream filtering, so we're going to use more conservative filters than typical
    gatk VariantFiltration \
        -R $refgenome \
        -V $rawindelsvcf \
        -O filtered_indels.vcf.gz \
        -filter-name "QD_filter" -filter "QD < 25.0" 
        
    gatk SelectVariants \
        --exclude-filtered \
        -V filtered_indels.vcf.gz \
        -O indels_filtered.vcf.gz
           
    gatk IndexFeatureFile \
        -I indels_filtered.vcf.gz
    """ 
}

process base_recal_boot_init{
    
    tag "$runAccession"
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict) 
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
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
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
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
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
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(bam), val(runAccession), path(table), path(snps), path(snp_idx), path(indels), path(indel_idx)


    output:
    tuple path('recal_3.bam'), val(runAccession)
    
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

process downstream_filter{
 
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    clusterOptions '--time 8:00:00 -A bharpur'
    publishDir "${params.savePath}/filtered_snps", mode: 'copy'

    input:
    tuple path(refgenome), path(refindex), path(dict)
    tuple path(combinedvcf), path(vcfindex)
    
    output:
    path('snps_filtered.vcf.gz')
    path('indels_filtered.vcf.gz')

    script:
    """
    ##to make the inputs
    gatk SelectVariants \
        -R $refgenome \
        -V $combinedvcf \
        -select-type SNP \
        -O raw_snps.vcf.gz
        
    gatk SelectVariants \
        -R $refgenome \
        -V $combinedvcf\
        -select-type INDEL \
        -O raw_indels.vcf.gz
    
    #filter SNPs using default params
    gatk VariantFiltration \
        -V raw_snps.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "SOR > 3.0" --filter-name "SOR3" \
        -filter "FS > 60.0" --filter-name "FS60" \
        -filter "MQ < 40.0" --filter-name "MQ40" \
        -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
        -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
        -O snps_filtered.vcf.gz

    #filter Indels using default params
    gatk VariantFiltration \
        -V raw_indels.vcf.gz \
        -filter "QD < 2.0" --filter-name "QD2" \
        -filter "QUAL < 30.0" --filter-name "QUAL30" \
        -filter "FS > 200.0" --filter-name "FS200" \
        -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
        -O indels_filtered.vcf.gz
        
    #output filtered:
    vcftools --gzvcf snps_filtered.vcf.gz --recode --remove-filtered-all --out snps_filtered.gz
    vcftools --gzvcf indels_filtered.vcf.gz --recode  --remove-filtered-all --out indels_filtered.gz
    
    """
}





/*
========================================================================================
========================================================================================
PIPELINE EXECUTION CONTROL
*/
workflow{
    
    index_genome(ch_refgenome)
    
    Channel.fromFilePairs(params.fastqPath+params.fastqPattern)
        | view() \
        | set {fastq_secured}

    alignment(index_genome.out,fastq_secured) \
        | set{align_done}
    
    check_duplicates(align_done) \
        | set {dupl_check_complete}

    remove_duplicates(dupl_check_complete) \
        | set {dupl_removed}

    if( params.knownSites != null ) {
        ch_knownsites = Channel.value(file(params.knownSites))
        base_recal1(index_genome.out, ch_knownsites, dupl_removed) \
            | set {base_recal1}
        base_recal2(index_genome.out, ch_knownsites, base_recal1) \
            | set {base_recal2}
        base_recal3(index_genome.out, ch_knownsites, base_recal2) \
            | set {base_recal_done}
        haplotype_caller(index_genome.out, base_recal_done) 
    } else {
        haplotype_caller(index_genome.out, dupl_removed)        
    }
    
    haplotype_caller.out.gvcf.collectFile(name: 'vcfs.list', newLine: true) \
        | set {vcfslist}
        
    combine_gvcfs(index_genome.out, vcfslist)
    genotype_gvcfs(index_genome.out, combine_gvcfs.out) \
    | set {combined_vcf}
    
    // If no known-sites file is provided, we apply a harsh variant filter and use those high-confidence sites to perform base quality recalibration
    if( params.knownSites == null ) {
        
        select_snps(index_genome.out, genotype_gvcfs.out)
        filter_snps(index_genome.out, select_snps.out)

        select_indels(index_genome.out, genotype_gvcfs.out)
        filter_indels(index_genome.out, select_indels.out)
            
        // Note that we use the first() operator here to ensure that the snp/indel filter outputs are treated as values instead of queues
        base_recal_boot_init(
            index_genome.out, 
            dupl_removed,
            filter_snps.out.first(),
            filter_indels.out.first())
            
        base_recal_boot_1(index_genome.out, base_recal_boot_init.out)
        base_recal_boot_2(index_genome.out, base_recal_boot_1.out)
        base_recal_boot_3(index_genome.out, base_recal_boot_2.out)
        base_recal_boot_3.out.view()
        
        haplotype_caller_2(index_genome.out, base_recal_boot_3.out) 
        haplotype_caller_2.out.gvcf.collectFile(name: 'vcfs.list', newLine: true) \
        | set {vcfslist_2}
        
        combine_gvcfs_2(index_genome.out, vcfslist_2)
        genotype_gvcfs_2(index_genome.out, combine_gvcfs_2.out) \
        | set {combined_vcf}
        
        genotype_gvcfs_2.out.view()
    }
    
    // Downstream quality filtering of SNPs
    downstream_filter(index_genome.out, combined_vcf)

}