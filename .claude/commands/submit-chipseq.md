# Submit ChIP-seq Analysis Job

Submit a ChIP-seq analysis job to the RunAI cluster.

## Quick Command

```bash
runai workspace submit chipseq-<SAMPLE> \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=<DATA_PATH>,mount=/data,mount-propagation=HostToContainer \
  --host-path path=<OUTPUT_PATH>,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=<REFERENCE_PATH>,mount=/references,mount-propagation=HostToContainer \
  --environment REFERENCE=/references/<INDEX_NAME> \
  --environment CPU_CORES=16 \
  --environment READ_TYPE=PE \
  --environment IP_SAMPLE=<IP_SAMPLE_NAME> \
  --environment INPUT_SAMPLE=<INPUT_SAMPLE_NAME> \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

## Parameters to Replace

| Placeholder | Description | Example |
|-------------|-------------|---------|
| `<SAMPLE>` | Job name suffix | `col0-h3k4me3` |
| `<DATA_PATH>` | Path to FASTQ files on cluster | `/hpi/hpi_dev/users/eberrigan/chipseq/data` |
| `<OUTPUT_PATH>` | Path for outputs on cluster | `/hpi/hpi_dev/users/eberrigan/chipseq/outputs` |
| `<REFERENCE_PATH>` | Path to reference files | `/hpi/hpi_dev/references/arabidopsis` |
| `<INDEX_NAME>` | Bowtie2 index name | `TAIR10` |
| `<IP_SAMPLE_NAME>` | IP sample file prefix | `IP_col0_H3K4me3` |
| `<INPUT_SAMPLE_NAME>` | Input control prefix | `input_col0` |

## Environment Variables

| Variable | Description | Options |
|----------|-------------|---------|
| `REFERENCE` | Bowtie2 index path | Path inside container |
| `CPU_CORES` | Parallel threads | `8`, `12`, `16` |
| `READ_TYPE` | Sequencing type | `PE` (paired-end), `SE` (single-end) |
| `IP_SAMPLE` | IP sample prefix | File prefix without `_R1.fq.gz` |
| `INPUT_SAMPLE` | Input sample prefix | File prefix without `_R1.fq.gz` |

## Example: Arabidopsis H3K4me3 ChIP-seq

```bash
runai workspace submit chipseq-col0-h3k4me3 \
  --project talmo-lab \
  --image ghcr.io/salk-harnessing-plants-initiative/chipseq-singlecell-pipeline:latest \
  --cpu-core-request 16 \
  --cpu-memory-request 64G \
  --host-path path=/hpi/hpi_dev/users/eberrigan/chipseq/data,mount=/data,mount-propagation=HostToContainer \
  --host-path path=/hpi/hpi_dev/users/eberrigan/chipseq/outputs,mount=/outputs,mount-propagation=HostToContainer,readwrite \
  --host-path path=/hpi/hpi_dev/references/arabidopsis,mount=/references,mount-propagation=HostToContainer \
  --environment REFERENCE=/references/TAIR10 \
  --environment CPU_CORES=16 \
  --environment READ_TYPE=PE \
  --environment IP_SAMPLE=IP_col0_H3K4me3 \
  --environment INPUT_SAMPLE=input_col0 \
  --command -- /scripts/chipseq/ChIProfiler_M00_Basic.sh
```

## Monitor Job

```bash
# Check status
runai workspace list | grep chipseq

# View logs
runai workspace logs chipseq-col0-h3k4me3 -p talmo-lab --follow

# Detailed status
runai workspace describe chipseq-col0-h3k4me3 -p talmo-lab
```

## Expected Outputs

After successful completion, `/outputs` will contain:
- `*_sorted_rmdup_uniq.bam` - Deduplicated alignments
- `*_TagDir/` - HOMER tag directories
- `*.bedGraph.gz` - Normalized coverage tracks
- `*.tdf` - IGV visualization files
- `1_res_read_mapping_stat.txt` - Mapping statistics
- `*_testPeaks*.bed` - Called peaks at various thresholds

## Pipeline Steps

1. **Trimming** - Adapter removal with trim_galore
2. **Alignment** - Map to reference with bowtie2
3. **BAM Processing** - Sort, deduplicate, filter
4. **Tag Directories** - Create HOMER tag directories
5. **Track Generation** - bedGraph and TDF files
6. **Statistics** - Read mapping statistics
7. **Peak Calling** - Call peaks with multiple thresholds

## Troubleshooting

### Job fails with "file not found"
- Verify `DATA_PATH` contains FASTQ files
- Check file naming matches `IP_SAMPLE` and `INPUT_SAMPLE` prefixes
- Ensure paired-end files have `_R1.fq.gz` and `_R2.fq.gz` suffixes

### Out of memory
Increase memory:
```bash
--cpu-memory-request 96G
```

### Reference index not found
- Verify bowtie2 index exists at `REFERENCE_PATH`
- Index should have files like `*.1.bt2`, `*.2.bt2`, etc.

## Related Commands

- `/monitor-jobs` - Monitor job status
- `/cleanup-jobs` - Clean up completed jobs
