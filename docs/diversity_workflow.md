# Background
We begin our diversity work using rarefied QIIME-formatted sequence and abundance table artifacts were created as follows:
1. The cutadapt-trimmed, DADA2 denoised representative sequences ([Mangan.raw_linked_required.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.raw_linked_required.repSeqs.qza)) had bat-host sequences removed, as described in the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md) document.
2. The remaining non-bat representative sequences ([Mangan.nonbatASVs.table.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.nonbatASVs.table.qza)) were filtered with the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script to retain only those ASVs that met two criteria:  
  **A.** The ASV was classified to the Arthropoda Phylum  
  **B.** The ASV contained at least Family-name information (i.e. Genus or Species name were retained, those missing Order name were discarded)  

  That filtering process resulted in generating a list of taxa ([taxfiltd_ASVs_NTCincluded.txt](https://github.com/devonorourke/mysosoup/blob/master/data/taxonomy/taxfiltd_ASVs_NTCincluded.txt)) that met those criteria, but required additional filtering considerations because a few negative control samples had some sequence data in addition to the guano samples.
  
3. The [contamination_investigations.md](https://github.com/devonorourke/mysosoup/blob/master/docs/contamination_investigations.md) document concludes with our assertion that the negative control samples are not indicative of pervasive contamination, and that ASVs detected among both guano and controls should be retained. As part of that process, the `taxfiltd_ASVs_NTCincluded.txt` list of taxonomically filtered ASVs generated from the `sequence_filtering.R` script was used to filter the original sequence ([Mangan.nonbatASVs.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.nonbatASVs.repSeqs.qza)) and table ([Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza)) artifacts to create ASV-filtered table ([Mangan.wNTCasvs-filt.table.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table.qza)) and sequence ([Mangan.wNTCasvs.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.wNTCasvs.repSeqs.qza)) artifacts.


> A rarefied table ([Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza)) was created by subsampling without replacement (using a sampling depth of 10,000 sequences) to perform diversity estimates as part of the `contamination_investigations.md` workflow.

> `$META` represents the full path to the QIIME-formatted [qiime_meta.tsv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/qiime_meta.tsv) file

```
qiime feature-table filter-samples \
  --i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
  --m-metadata-file $META \
  --p-where "SampleType='sample'" \
  --o-filtered-table Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza
```


4. These data were subsequently clustered using a 98.5% identity threshold.
```
qiime vsearch cluster-features-de-novo \
--i-sequences \
--i-table \
--p-perc-identity 0.985 \
--p-threads 4 \
--o-clustered-table \
--o-clustered-sequences

```


We start our bat diet analyses by discarding the negative control samples from the `Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza` table artifact, and drop any sequence variants unique to those samples:


In addition, we created a separate file for rarefied data containing only guano sample, but removed additional samples that contained only a single ASV. There were just 3 of the original 297 guano samples that were removed by this additional filter:
```
qiime feature-table filter-samples \
  --i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
  --m-metadata-file $META \
  --p-where "SampleType='sample'" \
  --p-min-features 2 \
  --o-filtered-table Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza
```

Finally, we applied those same two filters to the unrarefied dataset:
```
qiime feature-table filter-samples \
  --i-table Mangan.wNTCasvs-filt.table.qza \
  --m-metadata-file $META \
  --p-where "SampleType='sample'" \
  --p-min-features 2 \
  --o-filtered-table Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza
```

The [Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza) file served as input in the alpha diversity measures, while the [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza) artifact was used in beta diversity and supervised learning analyses. We used the [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza) file in the separate supervised learning  also to compare how model accuracy differed based on whether the input data was rarefied or not.


# Alpha diversity
We used a [Alpha-HillEstimates.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/Alpha-HillEstimates.R) script to produce the Hill Number estimates at each artificial roost site. The [Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza) file served as input for this analysis.

This script generated a plot of diversity estimates faceted by Hill Number, as well as a series of statistical summaries determining if diversity estimates varied by the Month or Site a sample was associated with. Multifactorial [ANOVA](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/anovas), [Kruskal-Wallis](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/kruskal), and [Dunn test](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/dunn) summaries are available for each of the three Hill Number diversity estimates.

The same R script also generated a few animations (which we could not display in a static manuscript). These are available at [this directory](https://github.com/devonorourke/mysosoup/tree/master/figures/gifs). We've embedded these animations in the [animation_analyses.md](https://github.com/devonorourke/mysosoup/blob/master/docs/animation_analyses.md) documentation and provided a brief commentary for each visualization.


# Beta diversity
Multiple phylogenetic (Unifrac weighted and unweighted) and non-phylogenetic (Dice-Sorensen, Bray-Curtis, Morisita-Horn) distances were estimated in the [betadiver_work.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/betadiver_work.R) script using the [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza) rarefied data table. The script resulted in three groups of output figures, tables, or summaries:
1. We ordinated each distance measure using PCoA to examine if patterns of collection Month and Location associated with sample ASV composition.
2. The Vegan function Adonis (multi factorial PERMANOVA) tested for Site and Month main effects using our five distance estimates. Summaries are available for each measure [in this directory](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/adonis).
3. We also measured whether the significant main effects observed in the Adonis tests were due to dispersion, thus we also used the Vegan PERMDISP (multi factorial PERMANOVA) function to test if within-group distances to group centroid differ across groups. Site and Month effects were tested separately, and a Tukey's HSD significance test was performed for each group within each distance measure. Summaries are [available in this directory](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/permdisp).

# Machine learning
The QIIME 2 [classify-samples](https://docs.qiime2.org/2019.1/plugins/available/sample-classifier/classify-samples/?highlight=classify%20samples) plugin offers a means to test a cross-validated supervised learning classifier - in this case, a Random Forest classifier. The goal is to take an input ASV table and a set of metadata variables (in our case, Site and/or Location when/where guano was collected) and train the classifier on a subset of these data. The resulting decision trees of the classifier are then applied to the other subset of data not used in training, and the accuracy of the classifier is assessed based upon how frequently the class variable prediction matches the truth. The individual ASV features used in training the classifier are individually removed and the subsequent accuracy of the model lacking that information can be evaluated to the original model to determine each features relative importance to the predictive power of the classifier. This key feature enables us to determine which ASVs are more or less important at distinguishing among the class variables like Site or Month. Thus, while ordinations can help understand if the overall composition of samples differ by these variables, a supervised learning approach helps identify the individual ASVs that are driving the largest differences among these variables.

We used the [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza) and [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza) artifacts as inputs in our supervised learning model because we were not clear whether randomly subsampling the data would have an impact on the overall model accuracy predictions. Note that subsampling these samples to a depth of 10,000 sequences resulted in discarding 7 of the 297 guano samples, and removed 84 ASVs from the original ASV table. However, those ASVs represent just 0.03% of the overall read abundance, thus it's highly unlikely that those ASVs discarded will significantly contribute to differences in our class variables of interest.

We executed scripts in QIIME 2 using the default settings with one exception: we increased the number of estimators from 100 to 1000 to potentially increase our model accuracy by increasing the number of trees in each forest. Models for Site, SiteMonth, CollectionMonth, and BatchType were performed as follows:

> `$TABLE` represents the [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza) rarefied dataset   
> `$METAFILE` represents the full path to the QIIME-formatted [qiime_meta.tsv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/qiime_meta.tsv) file


```
qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column CollectionMonth --output-dir learn-Month \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column Site --output-dir learn-Site \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column SiteMonth --output-dir learn-SiteMonth \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column BatchType --output-dir learn-BatchType \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000
```

The `feature_importance.qza` files were exported as `.tsv` files for use in subsequent R analyses. The `*accuracy.tsv` files [in this directory](https://github.com/devonorourke/mysosoup/tree/master/data/MachineLearn/rarefy) served as inputs to generating the heatmap plots made with the [machine_learn_heatmaps.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/machine_learn_heatmaps.R) script. The `fi*.tsv` files are used in the [machine_learn_analyses.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/machine_learn_analyses.R) script to generate several figures used in the analysis.


One of our first considerations was determining if we should be using rarefied or unrarefied data as input for the Random Forest classifier. We therefore applied the same set of scripts as noted above, but switched the input table to the [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza) unrarefied table.
> `$TABLE` represents the unrarefied [Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table_noNegSamps_noSingleASVs.qza) table  
> `$METAFILE` represents the full path to the QIIME-formatted [qiime_meta.tsv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/qiime_meta.tsv) file

```
qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column CollectionMonth --output-dir learn-Month \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column Site --output-dir learn-Site \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column SiteMonth --output-dir learn-SiteMonth \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000

qiime sample-classifier classify-samples \
  --i-table "$TABLE" --m-metadata-file "$METAFILE" \
  --m-metadata-column BatchType --output-dir learn-BatchType \
  --p-optimize-feature-selection --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000
```

Files were exported as noted above. The question of whether it is more appropriate to use rarefied data or unrarefied data was addressed in the initial portion of the [machine_learn_heatmaps.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/machine_learn_heatmaps.R) script. What quickly became clear was that rarefying the data increases the proportion of accuracy to the features (ASVs) that are most significant to the model, though the rank order of the important features tends not to vary substantially between rarefied or unrarefied datasets. In other words, if the point of the classifier is to identify the ASVs that best discriminate between Site or Month in which a sample is collected, whether or not you rarefy the data table isn't going to substantially change the top candidate ASVs. However, this comes with the caveat that _many more_ ASVs are considered as contributing to some fraction of variation in the unrarefied setting, and this may be because of differences in sequencing depth. Much like distance estimates needing some type of normalization to account for differences in abundances due to sampling/sequencing depth (we rarefy to address this, for instance), it's possible that the supervised learning estimates are more variable when per-sample sequencing depths aren't accounted for.

We plotted the cummulative sum of the feature importance (the relative change in model accuracy when it is removed from the classifier training set) in [this plot](https://github.com/devonorourke/mysosoup/blob/master/figures/feature_importance_Rarefied.vs.Nonrarefied.png) and were struck by how many more features are required to account for the top 90% model accuracy (red dotted line) in the unrarefied dataset for the **Month** and **SiteMonth** variables.

![SL_rare_nonrare_comp](https://github.com/devonorourke/mysosoup/blob/master/figures/feature_importance_Rarefied.vs.Nonrarefied.png)

We observed that rarefying data reduces the number of ASVs that contribute to the majority of model accuracy, thus, our analyses that went into plotting how these ASVs changed over Month or Site were selected from the rarefied supervised learning feature importance tables.

# Common ASVs
We examined what ASVs were identified in a high proportion of samples using the QIIME 2 [core-features](https://docs.qiime2.org/2019.4/plugins/available/feature-table/core-features/) plugin.

```
qiime feature-table core-features \
--i-table /mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.wNTCasvs-filt.rarefied-table_noNegSamps_noSingleASVs.qza \
--p-min-fraction 0.1 \
--p-steps 10 \
--o-visualization corefeatures.wNTCasvs.noNegs.noSingleASVs.qzv
```

These data are available in the [corefeatures.wNTCasvs.noNegs.noSingleASVs.qzv](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qzv/core_features/corefeatures.wNTCasvs.noNegs.noSingleASVs.qzv) visualization file. We used the [core_features.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/core_features.R) script to process these data and retained summaries of the information at [this directory](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/core_features).
