# Quantification of exon skipping





## SUPPA2 annotation

In script `src/suppa2_generate_events.sh`, we run SUPPA2 to get SE (skipped exon) events annotation. It creates a tab file. In addition, it calls `generate_CI_events.R` and `generate_CE_events.R` to extract constitutively included (CI) and excluded (CE) exons from the same genes as the SE. Finally, it makes two bed files from the tab file, to be used for quantification.

If needed, run `src/prep_alignments_bsn12.sh` (that in turns uses `src/combine_sj_files.R`) to assemble the technical replicates into a scratch directory. Only needed if it wasn't already done recently.

Then `src/get_psi_for_SE.sh` counts inclusion and exclusion reads with `bedtools` and `grep`, starting with SE annotation bed, and calls `src/assemble_psi.R` to create the result `assembled_psi.tsv`. Note: takes about 2.5 minutes/sample (8h total), and requires up to 200 GB memory for some samples.

Finally (not on cluster), can use `R/format_SE_psi.R` to load `assembled_psi.tsv`, Alec's sc-bulk integration thresholds for filtering, and export `events_coordinates.tsv` and `PSI_quantifications.tsv`.



# Versions

* 230531: WS281, bsn9 alignments
* 230907: updated CE
* 230912: updated CE
* 231109: WS289, bsn12 alignments
* 231214: temporary export without filtering
* 240308: quantifications (with proper filtering) of WS289, bsn12
* 240906: correct bug in CI that created 38 incorrect events (e.g. `CI_610` had `exon_start` before `intron_start`)
* 240910: correct additional bug in CI_8 (incorrect intron)



## Older approaches

### Schafer (2015) method

Initially tried with [Schafer](https://currentprotocols.onlinelibrary.wiley.com/doi/10.1002/0471142905.hg1116s87) method on WS281. It was successful on a test sample (ADLr173), but two problems:
* very long SJ (>20 kb) that skip entire genes are counted as exclusion reads. This leads to some exonic_parts that appear with intermediate PSI even though they are actually always included. Problem alleviated partially by filtering out long junctions BEFORE counting them. But better (non-exclusive) way would be to ensure that the SJ has at least 1 end within the gene before counting it. I can't think of an easy way to do that in the CLI framework.
* these PSI include all types of events, very few exon skipping but quite a bit of alt 3/5' ends alt first/last exons etc.

So, not ideal, we would like to restrict ourselves to proper exon skipping.

Scripts: `src/makeSchaferIndex.sh`, `src/schafer_quantif.sh`

### MISO

As Tsplice paper did, we can try MISO with an exon-centric annotation (MISO can either run for fll transcripts or exon-centric, but in that case we need an annotation of the individual events).

MISO provides exon-centric [annotations for some organisms](https://miso.readthedocs.io/en/fastmiso/annotation.html), not worms. In that case they suggest a tool, [rnaseqlib](https://rnaseqlib.readthedocs.io/en/clip/#creating-custom-gff-annotations-for-miso).

Turns out rnaseqlib is old and unmaintained, installing and running it runs into problems, gave up. Instead, going with a more recent software, SUPPA2, that includes script for creating annotations. See top of the page.


### SUPPA2 miniconda

```
conda create --name SUPPA2
conda activate SUPPA2
conda install -c bioconda suppa
```

### bsn9 prep annotations

*Not necessary anymore for bsn12: the technical replicates are already merged, can use bams directly.*


First, if necessary, run `prep_alignments_bsn9.sh` to assemble the bams and sj files from bsn9, store them on scratch60. It uses `samtools` and `src/combine_sj_files.R`. Important: McCleary doesn't have access to /SAY anymore, so the first step is to transfer all the bam files to scratch60, then assemble them.

```
input_bams="/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_bams/"
input_sj="/SAY/standard/mh588-CC1100-MEDGEN/bulk_alignments/bsn9_junctions/"
out_dir="/home/aw853/scratch60/2022-11-03_bsn9/"
```

use it to merge BAMs from technical replicates (using `samtools`), and merge SJs from technical replicates (using `src/combine_sj_files.R`).

Note: this step takes ~ 10 min per sample, ~30 h total.
