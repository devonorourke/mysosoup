# Classifying sequences
We explored both global alignment and a Naive Bayes classifier strategy to assign taxonomic information to the ASVs retained in the dataset. Multiple alignment parameters were explored as explained below, while only default settings for the Bayesian classifier was used.

## Alignment approaches

Non-host sequences were initially classified with VSEARCH, requiring a 97% identity and 94% coverage:
> "$REF" refers to reference file path to taxonomy mapping file and sequences
> "$REPSEQ" refers to path to `Mangan_noBats.repSeqs.qza` file

```
qiime feature-classifier classify-consensus-vsearch \
--i-query "$REPSEQ" \
--i-reference-reads "$REF"/ref_seqs.qza \
--i-reference-taxonomy "$REF"/ref_taxonomy.qza \
--p-maxaccepts 0 \
--p-strand both \
--p-threads 24 \
--o-classification Mangan_tax_p97c94.qza
```

The resulting `.qza` artifact was exported and used in the `sequence_filtering.R` script for further processing.
```
qiime tools export --input-path Mangan_tax_p97c94.qza --output-path mangan_tax_p97c94
mv mangan_tax_p97c94/* .
mv taxonomy.tsv mangan_tax_p97c94.tsv
```



We also explored two additional alignment parameters that reduced either the percent identity requirement, or the query coverage value:
- one alignment reduced the percent identity from 0.97 to 0.94 while maintaining the query coverage at 0.94
- a second alignment reduced the percent identity to 0.92, and lowered the coverage to 0.9

However, we applied these filters only on those data that were classified as "Unassigned" from the output (above). The rationale was to attempt to see what proportion of Unclassified sequences were assigned as such because of the percent identity applied in the initial alignment strategy. Sequences that were, for example, aligned at 96% identity can be included in some downstream analyses when taxonomic rank is incorporated. To filter for those Unassigned sequences from the original classified dataset, we first created a file of the Unassigned ASV IDs, and then created a new representative sequence set:

```
echo '#OTUID' > Unassigned_ASVs_manganp97c94.txt
grep 'Unassigned' mangan_tax_p97c94.tsv | cut -f 1 >> Unassigned_ASVs_manganp97c94.txt

REPSEQ=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan_noBats.repSeqs.qza
qiime feature-table filter-seqs \
--i-data "$REPSEQ" \
--m-metadata-file Unassigned_ASVs_manganp97c94.txt \
--o-filtered-data Mangan_Unassigned_noBats.repSeqs.qza
```

These `Mangan_Unassigned_noBats.repSeqs.qza` sequences were then aligned using the same script noted above, substituting the `--p-perc-identity` (from 0.97 to 0.94 or 0.92) and `--p-query-cov` (to 0.9) parameters.


## Naive Bayes classifier

We also classified sequences using a Naive Bayes approach through the SciKit Learn package implemented in QIIME2. This required first training the database using the combined arthropod, chordate, and host reference sequences:

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads /mnt/lustre/macmaneslab/devon/guano/BOLDdb/boldCOI.all.plusHost.seqs.qza \
  --i-reference-taxonomy /mnt/lustre/macmaneslab/devon/guano/BOLDdb/boldCOI.all.plusHost.tax.qza \
  --o-classifier nbClassifer_allPlusHost.qza
```

The trained classifier was then used to classify the non-bat sequences using default parameters:
> "$REPSEQ" refers to the full file path for the representative non-bat sequences `Mangan_noBats.repSeqs.qza`
> "$CLSSFYR" refers to the full file path for the trained classifier `nbClassifer_allPlusHost.qza`

REPSEQ=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan_noBats.repSeqs.qza
CLSSFYR=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/classifier/nbClassifer_allPlusHost.qza

```
qiime feature-classifier classify-sklearn \
--i-reads "$REPSEQ" --i-classifier "$CLSSFYR" \
--p-n-jobs 1 --p-reads-per-batch 2000 \
--o-classification mangan.sktax.qza
```
