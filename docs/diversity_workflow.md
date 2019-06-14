# Background
We begin our diversity work using rarefied QIIME-formatted sequence and abundance table artifacts were created as follows:
1. The cutadapt-trimmed, DADA2 denoised representative sequences ([Mangan.raw_linked_required.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.raw_linked_required.repSeqs.qza)) had bat-host sequences removed, as described in the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md) document.
2. The remaining non-bat representative sequences ([Mangan.nonbatASVs.table.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza)) were filtered with the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script to retain only those ASVs that met two criteria:  
  **A.** The ASV was classified to the Arthropoda Phylum  
  **B.** The ASV contained at least Family-name information (i.e. Genus or Species name were retained, those missing Order name were discarded)  

  That filtering process resulted in generating a list of taxa ([taxfiltd_ASVs_NTCincluded.txt](https://github.com/devonorourke/mysosoup/blob/master/data/taxonomy/taxfiltd_ASVs_NTCincluded.txt)) that met those criteria, but required additional filtering considerations because a few negative control samples had some sequence data in addition to the guano samples.
3. The [contamination_investigations.md](https://github.com/devonorourke/mysosoup/blob/master/docs/contamination_investigations.md) document concludes with our assertion that the negative control samples are not indicative of pervasive contamination, and that ASVs detected among both guano and controls should be retained. As part of that process, the `taxfiltd_ASVs_NTCincluded.txt` list of taxonomically filtered ASVs generated from the `sequence_filtering.R` script was used to filter the original sequence ([Mangan.nonbatASVs.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.nonbatASVs.repSeqs.qza)) and table ([Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza)) artifacts to create ASV-filtered table ([Mangan.wNTCasvs-filt.table.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table.qza)) and sequence ([Mangan.wNTCasvs.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.wNTCasvs.repSeqs.qza)) artifacts.
4. A rarefied table ([Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza)) was created by subsampling without replacement (using a sampling depth of 10,000 sequences) to perform diversity estimates as part of the `contamination_investigations.md` workflow.

We start our bat diet analyses by discarding the negative control samples from the `Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza` table artifact, and drop any sequence variants unique to those samples:

> `$META` represents the full path to the QIIME-formatted [qiime_meta.tsv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/qiime_meta.tsv) file

```
qiime feature-table filter-samples \
  --i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
  --m-metadata-file $META \
  --p-where "SampleType='sample'" \
  --o-filtered-table Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza
```

In addition, we similarly dropped negative control samples and ASVs from the unrarefied dataset also. However, in this case we also dropped any sample that contained only a single ASV (as these confound distance estimates). There were just 3 of the original 297 guano samples that were removed by this additional filter:
```
qiime feature-table filter-samples \
  --i-table Mangan.wNTCasvs-filt.table.qza \
  --m-metadata-file $META \
  --p-where "SampleType='sample'" \
  --p-min-features 2 \
  --o-filtered-table Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza
```

The resulting [Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza]() and [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza]() tables were used in the subsequent diversity measures described below.

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
