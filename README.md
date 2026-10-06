# humann_nf

Minimal Nextflow DSL2 workflow for HUMAnN 3.9. It uses paired, host-filtered
FASTQs and existing MetaPhlAn profiles, runs samples in parallel, and produces
sample-level and cohort-level functional tables.

## Shared defaults

The shared HUMAnN reference and container paths are defined in `conf/base.config`.
The input samplesheet and output path must be supplied in a project-local
parameter file; they are intentionally not stored in this shared pipeline.

Samplesheet columns: `sample,fastq_1,fastq_2,taxonomic_profile`.
Reference paths and every HUMAnN option are defined in `conf/base.config`.
Do not edit `main.nf` for routine parameter changes.

## Run

```bash
cd /path/to/project

nextflow run /data/software/nextflow/humann_nf \
  -profile singularity \
  -params-file /path/to/project/humann.params.yml \
  -work-dir work_humann \
  -resume
```

For a run-specific configuration, copy `params.yml`, edit the copy, and pass
that file with `-params-file`. Infrastructure settings can also be overridden
with a separate Nextflow config using `-c run.config`.

## Parallelism

The defaults run at most 10 samples concurrently with 8 CPUs per sample. Change
`max_parallel_samples` and `humann_cpus` together so their product stays below
the CPUs available to Nextflow.

## Important parameters

- `prescreen_threshold`: MetaPhlAn abundance threshold used to select pangenomes.
- `nucleotide_*_threshold`: nucleotide alignment filtering.
- `translated_*_threshold`: DIAMOND translated-search filtering.
- `bypass_nucleotide_search`: skip ChocoPhlAn search; normally keep `false`.
- `bypass_translated_search`: skip UniRef90 search; normally keep `false`.
- `remove_stratified_output`: removes species contributions; keep `false` when
  species-specific functions are required.
- `humann_extra_args`: escape hatch for valid HUMAnN 3.9 arguments not otherwise
  represented as a named parameter.
- `input`: CSV samplesheet with absolute FASTQ and MetaPhlAn profile paths.

## Outputs

- `humann/<sample>/`: gene families, pathway abundance, pathway coverage, log
- `joined/`: cohort-wide raw tables
- `normalized/`: CPM or relative-abundance tables
- `split/`: separate stratified and unstratified tables
- `pipeline_info/`: Nextflow trace, report, timeline, and DAG

Species-specific contributions are retained in stratified gene-family and
pathway-abundance tables. These are species-level assignments, not direct
strain-level functional evidence.

`pathcoverage` is intentionally not renormalized because it is a 0-1 pathway
completeness score rather than an abundance measurement.

## Restart

Always use `-resume` with the same work directory. Nextflow then reruns only
failed or changed tasks.
