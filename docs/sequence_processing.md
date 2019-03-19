# Data import
We're going to import the raw fastq files into QIIME by creating a file that QIIME2 reads to import all individual files into a single zipped `.qza` archive file in QIIME format. We don't have a typical input file type, so we'll create a manifest file and import as paired-end data. See [their import tutorial](https://docs.qiime2.org/2018.11/tutorials/importing/#sequence-data-with-sequence-quality-information-i-e-fastq) for full details.

1. To create the manifest file:

Starting in whatever project directory available:
```
## start within the project's fastq file directory
## find and remove any file pairs that one or both _R1_ or _R2_ files with zero bytes of data
find . -type f -empty -delete
pwd > ../qiime/inputs/pwd.tmptxt
ls -1 | cut -f 1 -d '_' | uniq -d > ../qiime/inputs/keeplist.tmptxt
find . -type f | grep -f ../qiime/keeplist.tmptxt | sort -u | cut -d '/' -f 2 > ../qiime/inputs/filenames.tmptxt
cd ../qiime/inputs
## create sample names to paste (col1)
cut -f 1 -d "_" filenames.tmptxt > col1.tmptxt
wc -l col1.tmptxt | cut -f 1 -d ' ' > lines.tmptxt
## create directory col to paste (col2)
seq $(echo $(cat lines.tmptxt)) | xargs -Iz echo $(cat pwd.tmptxt) > dirpath.tmptxt
paste dirpath.tmptxt filenames.tmptxt -d "/" > col2.tmptxt
## create read direction col to paste (col3)
for i in $(cat filenames.tmptxt); do
if [[ $i == *_R1_* ]];
then
  echo "forward"
else
  echo "reverse"
fi;
done > col3.tmptxt
paste col1.tmptxt col2.tmptxt col3.tmptxt -d "," > manifest.tmptxt
head -1 col2.tmptxt | cut -f 8 -d '/' > libname.tmptxt
echo 'sample-id,absolute-filepath,direction' | cat - manifest.tmptxt > $(cat libname.tmptxt).manifest.file
rm *.tmptxt
```

2. Import the data:
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path Mangan.manifest.file \
  --output-path Mangan.demux.qza \
  --input-format PairedEndFastqManifestPhred33
```

# Primer and barcode trimming with Cutadapt
We'll be using the QIIME cutadapt plugin to remove the primer and barcode sequences from the raw reads. See their [documentation](https://docs.qiime2.org/2018.11/plugins/available/cutadapt/?highlight=cutadapt) - specific to our data structure, we're using the [trim-paired](https://docs.qiime2.org/2018.11/plugins/available/cutadapt/trim-paired/) option because our data is already demultiplexed.

```
qiime cutadapt trim-paired \
  --i-demultiplexed-sequences Mangan.demux.qza \
  --p-cores 24 \
  --p-front-f GGTCAACAAATCATAAAGATATTGG \
  --p-adapter-f GGATTTGGAAATTGATTAGTWCC \
  --p-front-r GGWACTAATCAATTTCCAAATCC \
  --p-adapter-r CCAATATCTTTATGATTTGTTGACC \
  --o-trimmed-sequences Mangan.trimd.qza
```

We can create a per-sample summary of the resulting read depths that were trimmed with cutadapt. This visualization also includes an interactive plot that illustrates the per-sample base quality scores that are used to inform the DADA2 filtering parameters; note we're subsampling from a maximum of 100,000 reads per sample to generate base quality information.

```
qiime demux summarize \
  --i-data Mangan.trimd.qza \
  --p-n 100000 \
  --o-visualization Mangan.trimd.qzv
```  

The resulting visualization can be inspected to demonstrate that our forward and reverse sequencing quality appear as expected: reads are generally of high quality through the first 180 bp (our COI amplicon itself), but then start to deteriorate in quality thereafter. Forward reads drop significantly around 230 bp, while reverse reads drop off earlier around 200 bp. The DADA2 [tutorial](https://benjjneb.github.io/dada2/tutorial.html) suggests we trim the last 10 bp of the reads where quality deteriorates because it helps improve the algorithms sensitivity to rare sequence variants. We'll trim at the full length of the amplicon, a few base pairs where the quality scores _start_ to drop: 175 bp. There is more than sufficient overlap to join paired ends following DADA2 denoising of single end data.

## Denosing with DADA2
Previous testing indicated that DADA2 effectively reduces the number of unique sequence variants (ASVs) that are produced by the inherent errors in sequence base assignment generated during an Illumina run; these rare errors ultimately inflate the number of distinct amplicons in the dataset. DADA2 discards singleton sequences by default, but retains doubletons, etc.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs Mangan.trimd.qza \
  --p-trunc-len-f 175 \
  --p-trunc-len-r 175 \
  --p-n-threads 24 \
  --o-table Mangan.raw.table.qza \
  --o-representative-sequences Mangan.raw.repSeqs.qza \
  --o-denoising-stats Mangan.DADA2.denoisingStats.qza

qiime metadata tabulate \
  --m-input-file Mangan.DADA2.denoisingStats.qza \
  --o-visualization Mangan.DADA2.denoisingStats.qzv  
```

The resulting data following DADA2 denoising includes a set of representative sequences (`Mangan.raw.repSeqs.qza`), a matrix of the abundances of reads per ASV/Sample (`Mangan.raw.table.qza`), and summary files that illustrate the number of chimeric sequences that were discarded (`Mangan.DADA2.denoisingStats.qz*`). Data are further filtered to (1) remove host sequences, and (2) remove any remaining chimeras through reference-based detection.

# Filtering bat reads
We attempted to separate ASVs that contain host (bat) DNA from ASVs likely derived from diet by constructing a QIIME database containing a series of possible host reference sequences. These references included the focal species of _Myotis sodalis_ as well as a range of bats found across the Midwest and Northeast. In addition, we included other bat and bird species that had guano samples processed in our lab in previous projects; inclusion of these host sequences was done as a measure of contaminant check and removal. 
Representative ASVs from the `Mangan.raw.repSeqs.qza` file were queried for alignment to host reference sequences:

```
qiime quality-control exclude-seqs \
  --i-query-sequences ../inputs/Mangan_ASVseq.qza \
  --i-reference-sequences ../refs/batref_seqs.qza \
  --p-perc-identity 0.85 \
  --o-sequence-hits Mangan_hits.qza \
  --o-sequence-misses Mangan_misses.qza \
  --p-method vsearch
```
