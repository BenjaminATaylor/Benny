# nf_bennyben

A Nextflow pipeline for trimming, aligning, deduplicating, and joint variant-calling short-read sequencing data against a reference genome. Originally written for *Apis mellifera* by Benny Goller (Purdue College of Agriculture Data Services, October 2022); substantially rewritten and extended for hymenopteran genomics by Benjamin A. Taylor (Harpur Lab, 2024–2026).

## What it does

Given paired-end FASTQ files and a reference genome, the pipeline produces:

- Per-sample QC reports (fastp, FastQC, MultiQC, qualimap)
- Sorted, deduplicated, read-group-tagged BAMs
- Per-sample gVCFs (HaplotypeCaller)
- A jointly-called, hard-filtered, recoded multi-sample VCF (SNPs and indels)

### Processes, in order

| # | Process | Tool(s) | Purpose |
|---|---|---|---|
| 1 | `index_genome` | samtools, GATK | Build `.fai` and `.dict` for the reference |
| 2 | `fastp` | fastp | Adapter and quality trimming |
| 3 | `fastqc` + `multiqc` | FastQC, MultiQC | Read-level QC reports |
| 4 | `alignment` | NextGenMap + samtools sort | Produce sorted BAM |
| 5 | `check_duplicates` | Picard AddOrReplaceReadGroups | Add read-group tags |
| 6 | `remove_duplicates` | GATK MarkDuplicates | Remove sequencing duplicates |
| 7 | `qualimap` + `qualimap_collate` | qualimap | Alignment QC + collated report |
| 8 | *(optional)* `base_recal1/2/3` | GATK BQSR | Base-quality recalibration; only runs if `--knownSites` is provided |
| 9 | `haplotype_caller` | GATK HaplotypeCaller (gVCF mode) | Per-sample gVCF |
| 10 | `genomicsdb_import` | GATK GenomicsDBImport | Per-interval sample combine (parallelized) |
| 11 | `genotype_gvcfs` | GATK GenotypeGVCFs | Per-interval joint genotyping |
| 12 | `gather_vcfs` | GATK GatherVcfs | Recombine intervals in reference `.dict` order |
| 13 | `downstream_filter` | GATK SelectVariants + VariantFiltration | Split SNPs/indels and apply hard filters |
| 14 | `recode_vcfs` | vcftools | Mask genotypes with `DP < 3`; preserve INFO column |
| 15 | `fill_tags` | bcftools `+fill-tags` | Refresh AC/AN/AF/NS/F_MISSING; output bgzip + tabix |

## Parameters

### Required

- `--fqPattern <glob>` — Absolute path glob for paired FASTQs, e.g. `'/path/to/data/*_{1,2}.fq.gz'`. The literal `{1,2}` defines the pair.
- `--savePath <dir>` — Where published outputs are written. Subdirectories are created automatically (`raw_bams/`, `final_bam_files/`, `fastqc/`, `multiqc/`, `qualimap/`, `raw_snps/`, `filtered_snps/`).
- `--refGenome <file.fna>` — Absolute path to the reference FASTA.

### Optional

- `--intervals <file>` — Newline-separated list of intervals (one per line) for parallel joint genotyping. If omitted, intervals are derived from each contig in the reference `.dict`.
- `--knownSites <file.vcf.gz>` — A VCF of known variants. If supplied, BQSR (three iterations) runs before HaplotypeCaller. If omitted, BQSR is skipped.

## Requirements

- Nextflow 25.x or newer
- Singularity (or Apptainer) — every process runs in a container; no local tool installs needed
- A SLURM cluster (the default executor; configurable in `nextflow.config`)

## Running on a SLURM cluster

The repo ships a minimal `nextflow.config`. Cluster-specific options (account, partition, Tower token) should go in a personal config that you point Nextflow at with `-c`. Example for Purdue Negishi (`~/config/nextflow.config`):

```groovy
singularity.enabled = true
tower.accessToken = 'YOUR_TOWER_TOKEN'    // optional, only needed for -with-tower

executor {
    name = 'slurm'
    queueSize = 10
}

process {
    cache = 'lenient'
    clusterOptions = '-A bharpur -p cpu'  // Negishi requires explicit -A and -p
    errorStrategy = 'retry'
}
```

Then invoke:

```bash
nextflow \
  -c ~/config/nextflow.config \
  run /depot/bharpur/apps/nf_bennyben/main.nf \
  --fqPattern '/depot/bharpur/data/projects/example/data/*_{1,2}.fq.gz' \
  --savePath /depot/bharpur/data/projects/example/output \
  --refGenome /depot/bharpur/data/ref_genomes/VMAN/GCF_014083535.2_V.mandarinia_Nanaimo_p1.0_genomic.fna \
  -w /scratch/$USER/example_work \
  -bg -resume -with-tower
```

### Useful Nextflow flags

- `-bg` — run in the background; safe to close the terminal.
- `-resume` — pick up from cached tasks. Required for any incremental restart.
- `-with-tower` — stream progress to the Seqera platform. Requires `TOWER_ACCESS_TOKEN` in the environment or `tower.accessToken` in a config.
- `-w <dir>` — set the `work/` directory. **Recommended**: point at scratch (e.g. `/scratch/$USER/...`) so heavy intermediates don't accumulate on depot. Do not change between runs you want to `-resume` — the cache lives there.

## Troubleshooting

- **Logs**: `.nextflow.log` in the launch directory records every task. Search for `ERROR` or `terminated with an error` to find the failure.
- **Per-task logs**: a failure message points at a `work/XX/YYYYYY...` directory. Inside, the most useful files are:
  - `.command.log` — combined stdout/stderr from the task
  - `.command.run` — the wrapper SLURM script
  - `.exitcode` — the task's exit code
  - `.command.sh` — the actual command Nextflow generated
- **Session lock errors on `-resume`**: a previous run was killed without cleanup. Confirm nothing else is holding it (`lsof <work>/.nextflow/cache/*/db/LOCK` should be empty), then `rm` the LOCK file and rerun.
- **Stopping a running pipeline**:
  - Foreground: Ctrl-C.
  - Background: `kill $(cat .nextflow.pid)`, then `scancel` any in-flight SLURM jobs.
- **SLURM jobs visible**: `squeue -u $USER` lists your in-flight per-task jobs (each is tagged with its process and sample/interval).
