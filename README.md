<h1>Readme</h1>

Created by B. Goller, Purdue College of Agriculture Data Services, October 2022

Updated by B. A. Taylor, Harpur Lab, May 2024 (compare main.nf to main_nodownload.nf)

<h2>Overview</h2>
This Nextflow pipeline is designed to automate downloading and processing of Honey bee (Apis mellifera) sequence data. Taking an input list of desired sequences (local or NCBI/EBI), the pipeline then aligns those sequence to a reference genome (honey bee by default), does some quality control, and produces a final variant call format file for each sequence.

The Nextflow pipeline contains the following processes (in sequential order):
1)	fastq_setup: check whether one or two fastq files need to be downloaded
2)	prepFastq: prepare the one or two links for downloading
3)	download_fastq_links: download the fastq files needed
4)	alignment: NextGenMap + samtools sort to produce BAM file from single or paired end reads downloaded previously
5)	check_duplicates: gatk/Picard AddOrReplaceReadGroups
6)	remove_duplicates: gatk MarkDuplicates with –REMOVE_SEQUENCING_DUPLICATES set to ‘true’
7)	base_recal1/2/3: three iterations of gatk ApplyBQSR, gatk BaseRecalibrator, and gatk AnalyzeCovariates. Read more about this process here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR. The basic idea is that we identify sources of systematic bias in the rate at which mismatched calls (against the reference) are found, (for example, do we tend to find an overabundance of variants at the ends of reads), and then we adjust the base quality scores accordingly. If we know where we *should* be finding variants, we can supply a vcf of known sites to be masked from this process. In this pipeline, by default we use a VCF of known AMEL SNPs. But what if we don't have a set of known variants? In that case, we might wish to perform some kind of bootstrapping, e.g. first genotyping, then running the original data through the BSQR process using the most high-confidence SNPs from the first pass as known sites. This can be done multiple times until the data seem to converge. That capability isn't currently implemented, but could be!
2024 Update by Ben Taylor: bootstrapping is now implemented!
8)	haplotype_caller: gatk HaplotypeCaller to generate VCF, using default parameters. 
9)	alignment_cleanup: delete fastq and other intermediate files (bam_files, updated_bam_files, and recal_bam_files)

<h2>Required parameters</h2>
--fqPattern: An absolute filepath indicating the patterns of one or more pairs of fastq files. For example, if my raw data are called something like '/depot/popgenomes/USDAHornets/AGH432_2.fq.gz', then this parameter should be '/depot/popgenomes/USDAHornets/diploid_symlinks/*_{1,2}.fq.gz'
--savePath: An absolute filepath giving the directory in which the pipeline products should be saved
--refGenome: An absolute filepath giving the location of a reference sequence (usually an .fna file)

<h2>How to run on Purdue Bell</h2>

1)	There is limited space on your home directory on Bell (25GB) which will very quickly be filled by sequencing files and Nextflow logging – use your scratch directory (scratch/bell/userName or $RCAC_SCRATCH) where possible for temporary files
2)	Similarly, cd into your scratch directory (cd $RCAC_SCRATCH) before starting a Nextflow run so that pipeline logging can be stored (it will almost certainly be >25GB)
3)	Locations of Nextflow, and the pipeline components:
To get to the Nextflow program: /depot/bharpur/apps/nextflow
The bee pipeline is located at: /depot/bharpur/apps/nf_align_and_call_haplotype/main.nf
The config file for the pipeline: /depot/bharpur/apps/nf_align_and_call_haplotype/nextflow.config

<h3>Example Nextflow pipeline call</h3>

/depot/bharpur/apps/nextflow  \
-c /depot/bharpur/apps/nf_bennyben/nextflow.config \
run /depot/bharpur/apps/nf_bennyben/main.nf \
--fqPattern '/depot/bharpur/data/popgenomes/USDAHornets/diploid_symlinks/*_{1,2}.fq.gz' \
--savePath /depot/bharpur/data/popgenomes/USDAHornets/benny_full/output \
--refGenome /home/tayl1080/bharpur/data/ref_genomes/VMAN/GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_genomic.fna \
-bg -with-tower -resume

Breakdown:

/depot/bharpur/apps/nextflow runs the Nextflow app

-c /depot/bharpur/apps/nf_align_and_call_haplotype/nextflow.config tells Nextflow where the nextflow.config file is (sets some basic rules for how Nextflow operates)

run is the Nextflow command to run the pipeline

--fqPattern, --savePath, --refGenome: these all specify the data or output locations to use for this pipeline instance

[OPTIONAL] -bg runs Nextflow in the background so you can close your terminal without stopping the Nextflow pipeline (recommend this because the pipeline will take a long time to complete)

[OPTIONAL] -resume if you need to restart the pipeline (pick up where it stopped for instance) you can use the resume command. This will be ignored if you are starting a new run of the pipeline

[OPTIONAL] -with-tower lets you monitor the pipeline progress from a useful web-based UI. You need to have set this up beforehand: https://training.nextflow.io/latest/basic_training/seqera_platform/

<h2>Troubleshooting</h2>
The directory where you launch the Nextflow run matters…your logging files and /work/ directory that keeps track of your progress will both be saved there. Nextflow will run as a process (won’t show up in the jobs queue) and the pipeline will automatically submit individual jobs for each sequence (tagged with the accession number) for each part of the pipeline. 
Look at the .nextflow.log files to figure out where in the pipeline you are. If you need to dive deeper into a particular job, the log files will tell you what subdirectory in /work/ you should look at – each job will be recorded in a unique subdirectory within /work/. Then you can look at the .command.log, .exitcode, .command.run, etc files within that unique subdirectory to figure out what is happening/happened during the job.
See what jobs are currently running on the bharpur node: squeue -A bharpur (or for your username: squeue -u userName)

How to stop Nextflow once it starts running:
If you are not using the -bg (background) option: Nextflow is easy to stop, just close the Terminal or Cntrl/Command – C

When running with the -bg (background) it is harder to stop a Nextflow pipeline once you start it…
You can cancel any individual job using: scancel jobID#
The actual pipeline can be canceled using: kill processID#
You can find the processID# by looking at the .nextflow.pid file in the directory where you launched the Nextflow job (where the /work/ directory is)
