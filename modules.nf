/*
========================================================================================
    HAPLOTYPE CALLER
========================================================================================
*/

process haplotype_caller{
    
    time '1d'
    memory '10 GB'
    tag "$runAccession"
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    
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

process genomicsdb_import {

    time '12h'
    memory { 32.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    maxRetries 2
    tag "$interval"
    container "docker.io/broadinstitute/gatk:4.5.0.0"

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple path(vcfslist), val(interval)

    output:
    tuple val(interval), path("genomicsdb_workspace")

    script:
    """
    V_FLAGS=\$(awk '{printf "-V %s ", \$1}' ${vcfslist})
    gatk GenomicsDBImport \\
        \${V_FLAGS} \\
        --genomicsdb-workspace-path genomicsdb_workspace \\
        -L ${interval} \\
        --reader-threads 1 \\
        --batch-size 50
    """
}

process genotype_gvcfs {

    time '12h'
    memory { 32.GB * task.attempt }
    errorStrategy { task.exitStatus == 137 ? 'retry' : 'terminate' }
    maxRetries 2
    tag "$interval"
    container "docker.io/broadinstitute/gatk:4.5.0.0"

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    tuple val(interval), path(genomicsdb)

    output:
    tuple val(interval), path("genotyped_${interval}.vcf.gz"), path("genotyped_${interval}.vcf.gz.tbi")

    script:
    """
    gatk GenotypeGVCFs \\
        -R ${refgenome} \\
        -V gendb://${genomicsdb} \\
        -O genotyped_${interval}.vcf.gz

    gatk IndexFeatureFile \\
        -I genotyped_${interval}.vcf.gz
    """
}

process gather_vcfs {

    time '4h'
    memory '16 GB'
    container "docker.io/broadinstitute/gatk:4.5.0.0"
    publishDir "${params.savePath}/raw_snps", mode: 'copy'

    input:
    tuple path(refgenome), path(refindex), path(refdict)
    path vcfs
    path tbis

    output:
    tuple path('combined.vcf.gz'), path('combined.vcf.gz.tbi')

    script:
    """
    # Order inputs by reference order from the .dict (alphabetical can be wrong
    # when contig naming mixes prefixes, e.g. NC_/NW_ in NCBI assemblies where
    # NW scaffolds may precede the NC chromosome in the FASTA).
    # Assumes interval names match contig names (the default when intervals are
    # derived from the dict). Custom interval names that don't match a contig
    # would be skipped.
    I_FLAGS=""
    while IFS= read -r ctg; do
        f="genotyped_\${ctg}.vcf.gz"
        [ -f "\$f" ] && I_FLAGS="\${I_FLAGS} -I \$f"
    done < <(awk '/^@SQ/ { for (i=1;i<=NF;i++) if (\$i ~ /^SN:/) { sub("^SN:","",\$i); print \$i } }' ${refdict})

    gatk GatherVcfs \\
        \${I_FLAGS} \\
        -O combined.vcf.gz

    gatk IndexFeatureFile \\
        -I combined.vcf.gz
    """
}