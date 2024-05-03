Nextflow readme
Last edited by B. Goller, Purdue College of Agriculture Data Services on October 03, 2022

Overview:
This Nextflow pipeline is designed to automate downloading and processing of Honey bee (Apis mellifera) sequence data. Taking an input list of desired sequences (NCBI/EBI), the pipeline then aligns those sequence to a honey bee reference genome, does some quality control, and produces a final variant call format file for each sequence.

The Nextflow pipeline contains the following processes (in sequential order):
1)	fastq_setup: check whether one or two fastq files need to be downloaded
2)	prepFastq: prepare the one or two links for downloading
3)	download_fastq_links: download the fastq files needed
4)	alignment: NextGenMap + samtools sort to produce BAM file from single or paired end reads downloaded previously
5)	check_duplicates: gatk/Picard AddOrReplaceReadGroups
6)	remove_duplicates: gatk MarkDuplicates with –REMOVE_SEQUENCING_DUPLICATES set to ‘true’
7)	base_recal1/2/3: three iterations of gatk ApplyBQSR, gatk BaseRecalibrator, and gatk AnalyzeCovariates
8)	haplotype_caller: gatk HaplotypeCaller to generate VCF
9)	alignment_cleanup: delete fastq and other intermediate files (bam_files, updated_bam_files, and recal_bam_files)

Inputs:
--database: the information for the input sequences for the pipeline (format is a .CSV database). Subsetting the included ‘genomics_dataset.csv’ is recommended.
--saveTemp: indicate the directory where the downloaded fastq files and intermediate products should be saved (the fastq files and intermediate products are DELETED at the end of the workflow)
--savePath: indicate the directory where the pipeline products should be saved

Other notes:
1.	Reference genome: the pipeline is set up to work on Apis mellifera only – the reference genome has not been included as an argument and is therefore specified in the code directly. A future version could make the reference genome flexible.
2.	Similarly, a check for Apis mellifera sequences only is also included in the workflow. All database entries that do not have the ‘7460’ tax_id code are filtered out and ignored.

How to run on Purdue Bell:
Considerations:
1)	There is limited space on your home directory on Bell (25GB) which will very quickly be filled by sequencing files and Nextflow logging – use your scratch directory (scratch/bell/userName or $RCAC_SCRATCH) where possible for temporary files
2)	Similarly, cd into your scratch directory (cd $RCAC_SCRATCH) before starting a Nextflow run so that pipeline logging can be stored (it will almost certainly be >25GB)
3)	Make sure your “database” file with the fastq information is correct (contents and you are correctly pointing to the file) otherwise you will spend lots of time processing files you didn’t really want
4)	Locations of Nextflow, and the pipeline components:
To get to the Nextflow program: /depot/bharpur/apps/nextflow
The bee pipeline is located at: /depot/bharpur/apps/nf_align_and_call_haplotype/main.nf
The config file for the pipeline: /depot/bharpur/apps/nf_align_and_call_haplotype/nextflow.config

Example Nextflow pipeline call:
/depot/bharpur/apps/nextflow -C /depot/bharpur/apps/nf_align_and_call_haplotype/nextflow.config run /depot/bharpur/apps/nf_align_and_call_haplotype/main.nf --database /depot/bharpur/data/popgenomes/nextflow/PRJNA729035/genomics_dataset_PRJNA729035.csv --saveTemp $RCAC_SCRATCH --savePath /depot/bharpur/data/popgenomes/nextflow/PRJNA729035/ -bg -resume
Breakdown:
/depot/bharpur/apps/nextflow runs the Nextflow app
-C /depot/bharpur/apps/nf_align_and_call_haplotype/nextflow.config tells Nextflow where the nextflow.config file is (sets some basic rules for how Nextflow operates)
run is the Nextflow command to run the pipeline
/depot/bharpur/apps/nf_align_and_call_haplotype/main.nf is the Honey bee pipeline you want to run
--database /depot/bharpur/data/popgenomes/nextflow/PRJNA729035/genomics_dataset_PRJNA729035.csv points the pipeline to the CSV format database with study information (tells the pipeline what sequences to download, etc)
--saveTemp $RCAC_SCRATCH tells Nextflow to save the “temporary” files on your scratch space on Bell
--savePath /depot/bharpur/data/popgenomes/nextflow/PRJNA729035/ tells Nextflow to save the final output of the pipeline to a folder on the bharpur depot drive (in this case)
[OPTIONAL] -bg runs Nextflow in the background so you can close your terminal without stopping the Nextflow pipeline (recommend this because the pipeline will take a long time to complete)
[OPTIONAL] -resume if you need to restart the pipeline (pick up where it stopped for instance) you can use the resume command. This will be ignored if you are starting a new run of the pipeline

Troubleshooting:
The directory where you launch the Nextflow run matters…your logging files and /work/ directory that keeps track of your progress will both be saved there. Nextflow will run as a process (won’t show up in the jobs queue) and the pipeline will automatically submit individual jobs for each sequence (tagged with the accession number) for each part of the pipeline. 
Look at the .nextflow.log files to figure out where in the pipeline you are. If you need to dive deeper into a particular job, the log files will tell you what subdirectory in /work/ you should look at – each job will be recorded in a unique subdirectory within /work/. Then you can look at the .command.log, .exitcode, .command.run, etc files within that unique subdirectory to figure out what is happening/happened during the job.
See what jobs are currently running on the bharpur node: squeue -A bharpur (or for your username: squeue -u userName)

How to stop Nextflow once it starts running:
If you are not using the -bg (background) option: Nextflow is easy to stop, just close the Terminal or Cntrl/Command – C

When running with the -bg (background) it is harder to stop a Nextflow pipeline once you start it…
You can cancel any individual job using: scancel jobID#
The actual pipeline can be canceled using: kill processID#
You can find the processID# by looking at the .nextflow.pid file in the directory where you launched the Nextflow job (where the /work/ directory is)
