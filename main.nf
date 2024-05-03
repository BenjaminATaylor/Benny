#!/usr/bin/env nextflow
nextflow.enable.dsl=2
params.database = "/depot/bharpur/apps/nextflow_align_and_call_haplotype/genomics_dataset.csv"
//where to save the downloaded fastq files and initial versions of bam files (temporary, removed at the end)
params.saveTemp = "/depot/bharpur/data/popgenomes/nextflow"
//where to save the output files generated in the analyses
params.savePath = "/depot/bharpur/data/popgenomes/nextflow/"

//reference genome
params.refGenome = "/depot/bharpur/data/ref_genomes/AMEL/Amel_HAv3.1_genomic.fna"

//NCBI/EBI database taxonomic id for target species
//default is Apis mellifera
params.taxID = "7460"

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
//ORIGINAL VERSION DEVELOPED BY BENNY GOLLER, OCTOBER 2022
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
    SETUP FASTQ MANAGEMENT FOR ALIGNMENT
========================================================================================
*/

process fastq_setup{
    tag "$run_accession"
    module 'bioinfo:biopython/2.7.12'

    input:
    tuple val(run_accession), val(tax_id), val(scientific_name), val(fastq_ftp)

    output:
    stdout emit: var

    script:
    """
    #!/group/bioinfo/apps/apps/miniconda-4.10/bin/python
    #check to make sure we want to process this file

    fastq_all = '${fastq_ftp}'.split(';')

    if(len(fastq_all)==2):
        #want to run with 2 fastq files
        print(fastq_all[0])
        print(fastq_all[1])
        print('${run_accession}')

    else:
        #want to run with single fastq
        print(fastq_all[0])
        print("NULL")
        print('${run_accession}')
    """ 
}

/*
========================================================================================
    SETUP NEEDED DIRECTORIES, FASTQ DOWNLOAD AND ALIGNMENT
========================================================================================
*/

process download_fastq_links{
    tag "$runAccession"
    module 'bioinfo:samtools'
    clusterOptions '--ntasks 1 --time 8:00:00 -A bharpur'
    errorStrategy 'ignore'

    input:
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)

    script:
    if(fastq2 != "NULL")
        """
        mkdir -p ${params.saveTemp}/fastq_files
        mkdir -p ${params.saveTemp}/bam_files
        mkdir -p ${params.saveTemp}/updated_bam_files
        mkdir -p ${params.savePath}/final_bam_files
        mkdir -p ${params.savePath}/data_tables
        mkdir -p ${params.saveTemp}/recal_bam_files
        mkdir -p ${params.savePath}/recal_plots
        mkdir -p ${params.savePath}/raw_snps
        wget $fastq1 -O ${params.saveTemp}/fastq_files/"$runAccession"_1.fastq.gz
        wget $fastq2 -O ${params.saveTemp}/fastq_files/"$runAccession"_2.fastq.gz
        """
    else
        """
        mkdir -p ${params.saveTemp}/fastq_files
        mkdir -p ${params.saveTemp}/bam_files
        mkdir -p ${params.saveTemp}/updated_bam_files
        mkdir -p ${params.savePath}/final_bam_files
        mkdir -p ${params.savePath}/data_tables
        mkdir -p ${params.saveTemp}/recal_bam_files
        mkdir -p ${params.savePath}/recal_plots
        mkdir -p ${params.savePath}/raw_snps
        wget $fastq1 -O ${params.saveTemp}/fastq_files/"$runAccession".fastq.gz
        """
}

process alignment{
    tag "$runAccession"
    module 'bioinfo:samtools'
    clusterOptions '--ntasks 1 --time 1-0:00:00 --mem=4G -A bharpur'
    errorStrategy 'ignore'

    input:
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)

    script:
    if(fastq2 != "NULL")
        """
        /depot/bharpur/apps/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm -r ${params.refGenome} -1 ${params.saveTemp}/fastq_files/"$runAccession"_1.fastq.gz -2 ${params.saveTemp}/fastq_files/"$runAccession"_2.fastq.gz | samtools sort > ${params.saveTemp}/bam_files/"$runAccession".bam
        """
    else
        """
        /depot/bharpur/apps/NextGenMap-0.5.2/bin/ngm-0.5.2/ngm -r ${params.refGenome} -q ${params.saveTemp}/fastq_files/"$runAccession".fastq.gz | samtools sort > ${params.saveTemp}/bam_files/"$runAccession".bam
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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)

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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)


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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)


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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)


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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)


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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    tuple val(fastq1), val(fastq2), val(runAccession)


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
    tuple val(fastq1), val(fastq2), val(runAccession)

    output:
    val runAccession

    script:
    """
    rm ${params.saveTemp}/fastq_files/"$runAccession"*.fastq.gz
    rm ${params.saveTemp}/bam_files/"$runAccession".bam                        
    rm ${params.saveTemp}/updated_bam_files/"$runAccession"_updated.bam
    rm ${params.saveTemp}/recal_bam_files/"$runAccession"_recal_3*
    """
}
/*
========================================================================================
    PROCESS FASTQ LINKS FOR DOWNLOAD
========================================================================================
*/

process prepFastq{
    input:
    val pystr_out
    
    output:
    val fastq_split
    
    exec:
    fastq_split = pystr_out.split('\n')
    //println fastq_split
}

/*
========================================================================================
========================================================================================
PIPELINE EXECUTION CONTROL
*/
workflow{
    Channel.fromPath(params.database) \
        | splitCsv(header:true) \
        | map { row-> [row.run_accession, row.tax_id, row.scientific_name, row.fastq_ftp] } \
        | filter {list -> list[1] == params.taxID} \
        //| view() \
        | set {honeybee_records}
    
    fastq_setup {honeybee_records} \
        | set {fastq_to_align}

    prepFastq(fastq_to_align) \
        | map {val -> tuple(val[0], val[1], val[2])} \
        | view() \
        | set {fastq_for_download}

    download_fastq_links(fastq_for_download) \
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