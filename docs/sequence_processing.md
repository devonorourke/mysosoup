# Introduction
All information created for this project is available at [our Github repo](https://github.com/devonorourke/mysosoup). Please visit that page for more information regarding data tables, visualizations, and code used to complete this work.

# Molecular work  
Guano was provided as both single sample and pooled samples. Samples were processed either in 96-well plates or as single tube extractions using Qiagen's PowerSoil kits. Negative controls were used in both plates and single tubes to evaluate potential for cross-contamination during extraction.
Arthropod COI sequences were amplified using custom primers which targets a 180 bp region of cytochrome oxidase subunit-1. The resulting PCR products were normalized using Applied Biosystem's SequalPrep kit. The resulting normalized PCR products were pooled in fixed volumes and submitted for sequencing.
> These primers contain unique forward and reverse barcode/indices as well as full length Illumina primer index and flow-cell binding sequence. The result of a single PCR run generates ready-to-sequence amplicons.

# sequencing at NAU
The pooled library of COI amplicons were sequenced on an Illumina MiSeq machine with 300 bp PE sequencing using V2 chemistry set for 600 cycles at Northern Arizona University's sequencing center on September 14, 2018. **14,490,123 raw reads** were generated

# Installation specifications
We used a series of functions in QIIME2 (version 2019.1) to process the sequence data; a distinct Conda environment was created and followed specifications described in [QIIME2 documentation](https://docs.qiime2.org/2019.1/install/).

# Primer and barcode trimming with Cutadapt
While a very similar set of commands to remove non-COI sequence data (ie. barcodes, primers, etc.) exist within a QIIME environment there is a bug that exists in rare instances where a read is trimmed and results in a length of 0 - in this instance the program fails. The standalone Cutadapt version allows for an additional minimum length parameter which ensures that any reads with 0 (or whatever minimum length is used) are discarded in the filtered data - we set a minimum of 100 bases as the required length to discard any short fragments (the `-m 100` flag). In addition, we're employing the linked-read approach [see relevant Cutadapt section here](https://cutadapt.readthedocs.io/en/stable/guide.html#linked-adapters) that requires that the only reads which we retain are those in which the primer pairs are recognized in the read (for R1 this is the forward primer and the reverse complement of the reverse primer; for R2, it's the reverse primer and the reverse complement of the forward primer) - this is the `--trimmed-only` flag passed in the argument below. This requirement reduces the total number of retained reads, but it increases the likelihood that the read data we ultimately join (and retain for downstream analyses) are derived from instances in which both primers annealed at the expected areas of the COI fragment being amplified. In our experience, the linked read approach can reduce the number of ASVs by over 20% (compared to strictly passing the primers in their respective expected postions); in addition, the `--trimmed-only` flag retaining _only_ linked reads can further reduce spurious ASVs.

```
RAWDIR=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/rawfq

for SAMPLE in $(ls $RAWDIR | cut -f 1 -d '_' | sort -u); do
  cutadapt --cores=24 \
  -a GGTCAACAAATCATAAAGATATTGG...GGATTTGGAAATTGATTAGTWCC \
  -A GGWACTAATCAATTTCCAAATCC...CCAATATCTTTATGATTTGTTGACC \
  -m 100 --trimmed-only \
  -o "$SAMPLE"_1_trimd.fastq.gz -p "$SAMPLE"_2_trimd.fastq.gz \
  "$RAWDIR"/"$SAMPLE"_L001_R1_001.fastq.gz "$RAWDIR"/"$SAMPLE"_L001_R2_001.fastq.gz;
done
```

> `GGTCAACAAATCATAAAGATATTGG` represents the expected 5'-end primer for R1 data (forward primer)  
> `GGATTTGGAAATTGATTAGTWCC` represents the expected 3'-end primer for R1 data (reverse complement reverse primer)  
> `GGWACTAATCAATTTCCAAATCC` represents the expected 5'-end primer for R2 data (reverse primer)  
> `CCAATATCTTTATGATTTGTTGACC` represents the expected 3'-end primer for R2 data (reverse complement forward primer)

Note that the resulting trimmed reads remain unpaired. We next import these unpaired data into QIIME2 for DADA2 denoising and paired-end joining. Prior to denoising these data are evaluated for per-base quality to determine what fragment lengths are retained/trimmed by DADA2.

# Data import
Trimmed fastq files were imported into QIIME by creating a (manifest) file that QIIME2 reads to import all individual files into a single zipped `.qza` archive file in QIIME format. We don't have a typical input file type, so we'll create a manifest file and import as paired-end data. See [their import tutorial](https://docs.qiime2.org/2018.11/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq) for full details. Data was imported as follows:

```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../Mangan.manifest_linked_required.file \
  --output-path Mangan.trimd_linked_required.qza \
  --input-format PairedEndFastqManifestPhred33
```

> an example of the data structure of the `Mangan.manifest_linked_required.file` is as follows:
```
sample-id,absolute-filepath,direction
6212017EGA1,/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/trimd_linkReqd_fq/6212017EGA1_1_trimd.fastq.gz,forward
6212017EGA1,/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/trimd_linkReqd_fq/6212017EGA1_2_trimd.fastq.gz,reverse
6212017EGA2,/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/trimd_linkReqd_fq/6212017EGA2_1_trimd.fastq.gz,forward
6212017EGA2,/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/trimd_linkReqd_fq/6212017EGA2_2_trimd.fastq.gz,reverse
```

We can create a per-sample summary of the resulting read depths that were trimmed with cutadapt. This visualization also includes an interactive plot that illustrates the per-sample base quality scores that are used to inform the DADA2 filtering parameters; note we're subsampling from a maximum of 10,000 reads per sample to generate base quality information.

```
qiime demux summarize \
  --i-data Mangan.trimd_linked_required.qza \
  --p-n 10000 \
  --o-visualization Mangan.required_trimd_perBase_viz.qzv
```  

The resulting visualization can be inspected to demonstrate that our forward and reverse sequencing quality appear as expected: reads are generally of high quality through the first 181 bp (our target COI region of the amplicon), but then immediately deteriorate in quality thereafter. In fact, most reads do not contain more than the expected 181 bp (91% of reads are of the expected length); the remaining reads that exceed this length are of dramatically lower quality. Because our amplicon is expected to be 181 bp in length, reads exceeding this length likely represent samples that lacked primers and were left untrimmed. The DADA2 [tutorial](https://benjjneb.github.io/dada2/tutorial.html) suggests we trim the last 10 bp of the reads where quality deteriorates because it helps improve the algorithms sensitivity to rare sequence variants. We'll trim at the full length of the amplicon, a few base pairs where the quality scores _start_ to drop: 175 bp. There is more than sufficient overlap to join paired ends following DADA2 denoising of single end data. Furthermore, trimming around the 180 bp mark should sufficiently discard those instances in which a read lacked a primer and was untrimmed.

## Denosing with DADA2
Previous testing indicated that DADA2 effectively reduces the number of unique sequence variants (ASVs) that are produced by the inherent errors in sequence base assignment generated during an Illumina run; these rare errors ultimately inflate the number of distinct amplicons in the dataset. DADA2 discards singleton sequences by default, but retains doubletons, etc.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs Mangan.trimd_linked_required.qza \
  --p-trunc-len-f 175 \
  --p-trunc-len-r 175 \
  --p-n-threads 24 \
  --o-table Mangan.raw_linked_required.table.qza \
  --o-representative-sequences Mangan.raw_linked_required.repSeqs.qza \
  --o-denoising-stats Mangan.DADA2.denoisingStats_linked_required.qza

qiime metadata tabulate \
  --m-input-file Mangan.DADA2.denoisingStats_linked_required.qza \
  --o-visualization Mangan.DADA2.denoisingStats_linked_required.qzv
```

The resulting data following DADA2 denoising includes a set of representative sequences (`Mangan.raw_linked_required.repSeqs.qza`), a matrix of the abundances of reads per ASV/Sample (`Mangan.raw_linked_required.table.qza`), and summary files that illustrate the number of chimeric sequences that were discarded (`Mangan*.qz*`). Notably, four samples no longer contain any reads - all four were extraction blanks (negative controls). While 11 blanks were created, the remaining 7 samples contain substantial number of reads and likely reflect instances in which a piece of guano was accidentally transferred into a well. These will be further investigated for potential downstream effects on possible contamination.

We remove these four samples from the `.qza` artifact given that they have no data:
```
qiime feature-table filter-samples \
--i-table Mangan.raw_linked_required.table.qza \
--p-min-features 1 \
--o-filtered-table Mangan.samp-filt.table.qza
```

We also remove any ASVs that were associated _only_ with those four samples:
```
qiime feature-table filter-seqs \
--i-data Mangan.raw_linked_required.repSeqs.qza \
--i-table Mangan.samp-filt.table.qza \
--o-filtered-data Mangan.samp-filt.repSeqs.qza
```
  > It turns out this didn't remove any further samples

These data represent the input into the next phase of data processing: identifying (and separating) host sequences from non-host samples using a series of databases and classification approaches. See the `classify_sequences.md` document for the next step of the pipeline; as noted in that documentation, classification required the creation of custom databases - these are further explained in separate documents (see `host_database.md` and `database_construction.md`).
