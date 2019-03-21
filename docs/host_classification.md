# Classifying sequences
Host sequences were initially classified with VSEARCH, requiring a 98% identity and 94% coverage:
> "$REF" refers to reference file path to taxonomy mapping file and sequences
> "$REPSEQ" refers to path to `Mangan_noBats.repSeqs.qza` file

REF=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/host_dbs
REPSEQ=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan_BatsOnly.repSeqs.qza

```
qiime feature-classifier classify-consensus-vsearch \
--i-query "$REPSEQ" \
--i-reference-reads "$REF"/boldCOI.all.plusHost.seqs.qza \
--i-reference-taxonomy "$REF"/boldCOI.all.plusHost.tax.qza \
--p-perc-identity 0.98 \
--p-query-cov 0.94 \
--p-maxaccepts 0 \
--p-strand both \
--p-threads 24 \
--o-classification host_tax_p98c94.qza
```

The resulting `.qza` artifact was exported and used in the `sequence_filtering.R` script for further processing.
```
qiime tools export --input-path host_tax_p98c94.qza --output-path host_tax_p98c94
mv mangan_tax_p98c94/* .
mv taxonomy.tsv mangan_tax_p98c94.tsv
```

We also explored requiring a greater percent identity to see if sequences would avoid any ambiguity in host specificity by increasing identity to 99%:

```
qiime feature-classifier classify-consensus-vsearch \
--i-query "$REPSEQ" \
--i-reference-reads "$REF"/boldCOI.all.plusHost.seqs.qza \
--i-reference-taxonomy "$REF"/boldCOI.all.plusHost.tax.qza \
--p-perc-identity 0.99 \
--p-query-cov 0.94 \
--p-maxaccepts 0 \
--p-strand both \
--p-threads 24 \
--o-classification host_tax_p99c94.qza
```
