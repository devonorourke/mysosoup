# Background
The [contamination_investigations.md](https://github.com/devonorourke/mysosoup/blob/master/docs/contamination_investigations.md) document describes efforts used to justify retaining ASVs for our bat diet analyses for such instances where sequence variants were detected in both negative control samples and guano samples. This process began by taking the host-filtered representative sequence variants generated at the output of the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md) document. We used the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script to identify which ASVs were present in both guano and negative control samples. Notably, we restricted our analyses to those ASVs that met two criteria:
  - 1) The ASV was classified to the Arthropoda Phylum
  - 2) The ASV contained at least Family-name information (i.e. Genus or Species name were retained, those missing Order name were discarded)

Part of that contamination investiation applied diversity estimates, and these estimates required rarefying data. We begin our diversity work using those rarefied QIIME sequences and abundance tables ([Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza)) that contained ASVs detected in both guano and negative control samples, though the control samples are discarded as are any ASVs that were unique to those control samples.


## Alpha diversity
We used a [Alpha-HillEstimates.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/Alpha-HillEstimates.R) script to produce the Hill Number estimates at each artificial roost site. We generated the plot of diversity estimates faceted by Hill Number, and generated a series of statistical summaries determining if diversity estimates varied by the Month or Site a sample was associated with. Multifactorial [ANOVA](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/anovas), [Kruskal-Wallis](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/kruskal), and [Dunn test](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/dunn) summaries are available for each of the three Hill Number diversity estimates. 

## Beta diversity
See `betadiv_analyses.R`. Includes NMDS plots and Adonis estimates.

## Machine learning
Generating a .qza table and repSeq file that reflects:
  1. ASVs we want (Arthropod only, Family-rank required); ASVs only in control samples removed
  2. Samples we want (Samples used in beta diversity estimates; dropping 3 samples with just one ASV)

```
echo sampleid > samples2drop.txt
echo 7272017EGC1 >> samples2drop.txt
echo 7272017EGC2 >> samples2drop.txt
echo 7272017HBA6 >> samples2drop.txt

qiime feature-table filter-samples \
  --p-exclude-ids \
  --i-table Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza \
  --m-metadata-file samples2drop.txt \
  --o-filtered-table tmp.input.qza

qiime feature-table filter-features \
  --i-table tmp.input.qza \
  --p-min-samples 1 \
  --o-filtered-table Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza

rm tmp.input.qza  
```

## Machine learning

Running tests for the following metadata: Site, SiteMonth, CollectionMonth, BatchType

```
READS=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza
METAFILE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/meta/qiime_meta.tsv

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column CollectionMonth --output-dir learn-Month \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column Site --output-dir learn-Site \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column SiteMonth --output-dir learn-SiteMonth \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column BatchType --output-dir learn-BatchType \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000
```

The `feature_importance.qza` files were exported as `.tsv` files for use in subsequent R analyses:
```
qiime tools export --input-path ./{$DIR}/feature_importance.qza --output-path fi
```

Also rerunning the same Machine-leanring tests, but inputting a data set that isn't rarefied. Same ASVs, same samples:
```
echo sampleid > moresamples2drop.txt
echo 7272017EGC1 >> moresamples2drop.txt
echo 7272017EGC2 >> moresamples2drop.txt
echo 7272017HBA6 >> moresamples2drop.txt
echo ExtractionNTC1S1 >> moresamples2drop.txt
echo ExtractionNTC7S55 >> moresamples2drop.txt
echo ExtractionNTC11S115 >> moresamples2drop.txt
echo ExtractionNTC3S19 >> moresamples2drop.txt
echo ExtractionNTC8S64 >> moresamples2drop.txt
echo ExtractionNTC4S28 >> moresamples2drop.txt
echo ExtractionNTC2S10 >> moresamples2drop.txt
echo blankS39 >> moresamples2drop.txt

qiime feature-table filter-samples \
  --p-exclude-ids \
  --i-table Mangan.wNTCasvs-filt.table.qza \
  --m-metadata-file moresamples2drop.txt \
  --o-filtered-table tmp.input.qza

qiime feature-table filter-features \
  --i-table tmp.input.qza \
  --p-min-samples 1 \
  --o-filtered-table Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza

READS=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza
METAFILE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/meta/qiime_meta.tsv

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column CollectionMonth --output-dir learn-Month \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column Site --output-dir learn-Site \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column SiteMonth --output-dir learn-SiteMonth \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$READS" --m-metadata-file "$METAFILE" \
  --m-metadata-column BatchType --output-dir learn-BatchType \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000
```


## Common ASVs
We examined what ASVs were identified in a high proportion of samples:

```
qiime feature-table core-features \
--i-table /mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza \
--p-min-fraction 0.1 \
--p-steps 10 \
--o-visualization corefeatures.wNTCasvs.noNegs.noSingleASVs.qzv
```
