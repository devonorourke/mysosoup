# Background

Prior to beginning diversity analyses, we completed the following steps:

- Raw sequence reads were trimmed with cutadapt and denoised with DADA2 into representative sequences resulting in the [Mangan.raw_linked_required.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.raw_linked_required.repSeqs.qza)) object, as described in the [sequence_processing.md](https://github.com/devonorourke/mysosoup/blob/master/docs/sequence_processing.md) workflow.    

- We identified bat-host sequences, as described in the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md) document. The host database design was described in the [host_database.md](https://github.com/devonorourke/mysosoup/blob/master/docs/host_database.md) document.  

- The non-bat representative sequences were investigated for contamination as described in the [contamination_investigations.md](https://github.com/devonorourke/mysosoup/blob/master/docs/contamination_investigations.md) document - we failed to detect extensive contamination either during DNA extraction or PCR amplification. See also the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script related to this contamination investigation.  

# Sequence processing prior to diversity estimates
We begin this diversity workflow using the original dereplicated representative sequences. We next perform the following actions prior to applying diversity measures:  

1. Filter samples  
All pooled samples and negative control samples are discarded, as well as any sequence feature associated exclusively with those controls.  

> `$META` represents the full path to the QIIME-formatted [qiime_meta.tsv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/qiime_meta.tsv) file
```
## filter table to retain only guano sampled as single pellets
qiime feature-table filter-samples \
  --i-table Mangan.raw_linked_required.table.qza \
  --m-metadata-file qiime_meta.tsv \
  --p-where "SampleType='sample' AND BatchType='single'" \
  --p-min-features 1 \
  --o-filtered-table Mangan.dada2_singles_table.qza

## drop any ASVs exclusive to negative controls
qiime feature-table filter-seqs \
  --i-data Mangan.raw_linked_required.repSeqs.qza \
  --i-table Mangan.dada2_singles_table.qza \
  --o-filtered-data Mangan.dada2_singles_seqs.qza
```

2. Cluster samples  

Remaining representative sequences are clustered at 98.5% identity using `qiime vsearch cluster-features-de-novo`  
```
qiime vsearch cluster-features-de-novo \
  --i-table Mangan.dada2_singles_table.qza \
  --i-sequences Mangan.dada2_singles_seqs.qza \
  --p-perc-identity 0.985 \
  --o-clustered-table Mangan.clust_p985_table.qza \
  --o-clustered-sequences mv 
```

> Dropping the negative control and pooled samples, coupled with clustering reduces the number of unique sequence variants to 1,936 sequence representatives. Notably, there can still be many suprious diet components we may want to discard, including host DNA and non-arthropod COI sequences. In addition, the classification accuracy of some of these sequences may be insufficient for our diversity estimates, so we'll first classify all the sequences and then filter appropriately.

3. Classify samples  

Clustered sequences are classified using a hybrid approach: exact alignments with VSEARCH are prioritized first, and any sequences not classified with this method are passed to the naive Bayes approach with `qiime feature-classifier classify-hybrid-vsearch-sklearn`.  
First, for VSEARCH: 
```
qiime feature-classifier classify-consensus-vsearch --p-threads 28 \
--i-query Mangan.clust_p985_seqs.qza \
--i-reference-reads bigCOI.derep.seqs.qza \
--i-reference-taxonomy bigCOI.derep.tax.qza \
--p-perc-identity 1.0 --p-query-cov 0.94 \
--o-classification Mangan.clust_p985_taxa_VsearchOnly_p100_c94.qza

qiime tools export --input-path Mangan.clust_p985_taxa_VsearchOnly_p100_c94.qza --output-path tmpdir_VsearchOnly_p100c94
```

Next, for naive BAYES. For simplicity sake, I reclassified all sequences, but filtered the exact VSEARCH classifications in a subsequent R script.
```
qiime feature-classifier classify-sklearn --p-threads 20 --p-no-prefilter \
--i-reads Mangan.clust_p985_seqs.qza \
--i-classifier nbClassifer_bigDB_2020-2.qza \
--o-classification Mangan.clust_p985_taxa_sklearnOnly.qza

qiime tools export --input-path Mangan.clust_p985_taxa_sklearnOnly.qza --output-path tmpdir_sklearnOnly
```

4. Filter classified samples  

The outputs of these two `tmpdir_*` directories contain the text versions of the taxonomy labels for each OTU. They are used in an R script, [filtering_taxa_fromVsearchSKlearn.R](https://raw.githubusercontent.com/devonorourke/mysosoup/master/scripts/r_scripts/filtering_taxa_fromVsearchSKlearn.R), so that we obtain a final list of feature IDs to retain in a taxonomy-based filter that ensures:
- we retain all positive VSEARCH matches (for 100% identity and at least 94% query coverage) that are arthropods with at least family-rank information
- we retain the remaining positive matches for sklearn that are arthropods with at least family-rank information

This R script produces a pair of (text file) inputs listing the `Feature ID` values needed to filter the VSEARCH and naive Bayes classifier .qza files so that we retain only these particular features we want (that is, those Feature IDs initially and preferentially from VSEARCH, and then from naive Bayes, given the minimum amount of taxonomic information). These files are those used with the `--m-ids-to-keep-file` parameter in the commands below:

- Filter the VSEARCH file to retain just the Feature IDs required:
```
qiime rescript filter-taxa \
  --i-taxonomy Mangan.clust_p985_taxa_VsearchOnly_p100_c94.qza \
  --m-ids-to-keep-file filtd_taxlist_vsearch.txt \
  --o-filtered-taxonomy Mangan.clust_p985_taxa_Filtd_vsearch.qza
```

- Filter the naive Bayes classifier samples to retain just the Feature IDs required:
```
qiime rescript filter-taxa \
  --i-taxonomy Mangan.clust_p985_taxa_sklearnOnly.qza \
  --m-ids-to-keep-file filtd_taxlist_sklearn.txt \
  --o-filtered-taxonomy Mangan.clust_p985_taxa_Filtd_sklearn.qza
```

- Combine the two taxonomy files and used the COMBINED lists of vsearch/sklearn-filtered data to filter the clustered seqs object
```
qiime feature-table merge-taxa \
  --i-data Mangan.clust_p985_taxa_Filtd_vsearch.qza Mangan.clust_p985_taxa_Filtd_sklearn.qza \
  --o-merged-data Mangan.clust_p985_taxa_Filtd_All.qza

cut -f 1 -d ',' filtd_tax_dataframe_ALL.csv > filtd_taxlist_all.txt

qiime feature-table filter-seqs \
  --i-data Mangan.clust_p985_seqs.qza \
  --o-filtered-data Mangan.clust_p985_seqs_Filtd.qza \
  --m-metadata-file filtd_taxlist_all.txt
```

- Use the COMBINED lists of vsearch/sklearn-filtered data to filter the clustered table object
```
qiime feature-table filter-features \
  --i-table Mangan.clust_p985_table.qza \
  --o-filtered-table Mangan.clust_p985_table_Filtd.qza \
  --m-metadata-file filtd_taxlist_all.txt
```

5. Rarefy samples  

We first evaluate what the appropriate sampling depth would be by summarizing the feature table to estimate bounds of sampling depth options, then generate a rearefaction curve around that sampling depth boundary: 
```
qiime feature-table summarize --i-table Mangan.clust_p985_table_Filtd.qza --o-visualization Mangan.clust_p985_tableSummary.qzv

qiime diversity alpha-rarefaction \
  --i-table Mangan.clust_p985_table_Filtd.qza \
  --p-min-depth 2000 --p-max-depth 12000 \
  --p-metrics shannon observed_otus \
  --o-visualization Mangan.clust_p985_rarefactionCurve.qzv
```

We ultimately decided on a rarefying depth of 10,000 sequences per sample:
```
qiime feature-table rarefy \
  --i-table Mangan.clust_p985_table_Filtd.qza \
  --p-sampling-depth 10000 \
  --o-rarefied-table Mangan.clust_p985_table_Rarefyd.qza
```

The rarefied OTU table, [Mangan.clust_p985_table_Rarefyd.qza](https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_table_Rarefyd.qza), is the file that is used for alpha and beta diversity estimates described below. 


# Diversity estimates 

## Alpha diversity
Observed richness and Shannon's entropy were for each guano sample and calculated in Vegan using a custom R script, [alpha_analyses_singleOnly.R](https://raw.githubusercontent.com/devonorourke/mysosoup/master/scripts/r_scripts/alpha_analyses_singleOnly.R). This script generated the diversity estimates shown in Figure 2, as well as the supplementary table[all_alpha_pairwise_Wilcoxon.csv](https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/text_tables/wilcoxon_pairwise/all_alpha_pairwise_Wilcoxon.csv), the resport of a Wilcoxon rank sum test comparing pairwise differences in each alpha diversity metric for all possible Site + Month combinations.


## Beta diversity
While the alpha diversity estimates required only the rarefied OTU table, the beta diversity tests incorporated some phylogenetic-weighted distance measures that required providing a rooted tree. We built this in QIIME as follows prior to running the beta diversity tests:
```
## align the clustered seqs with mafft and build a rooted tree with fasttree
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences Mangan.clust_p985_seqs_Filtd.qza \
  --output-dir Mangan_mafft-fasttree-output
## export the aligned tree 
qiime tools export \
--input-path ./Mangan_mafft-fasttree-output/rooted_tree.qza \
--output-path newick_tree_dir
## renamed 'tree.nwk' as 'Mangan.clust_p985_rooted_tree.nwk'
mv tree.nwk Mangan.clust_p985_rooted_tree.nwk
```

We used the output of the newick tree, [Mangan.clust_p985_rooted_tree.nwk](https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/trees/Mangan.clust_p985_rooted_tree.nwk), in following R script (below).

Four distance metrics were applied to the rarefied OTU table using a custom R script, [beta_analyses_singleOnly.R](https://raw.githubusercontent.com/devonorourke/mysosoup/master/scripts/r_scripts/beta_analyses_singleOnly.R) to evaluate if community composition was associated with either collection Site or Month:
- non phylogenetic measures that were unweighted (Bray-Curtis) or weighted (Bray-Curtis) with respect to sequence abundance, or
- phylogenetic measures that were weighted (Weighted Unifrac) or unweighted (Unweighted Unifrac)

The R script resulted in three summary outputs:
1. The ordinations of each distance measure using PCoA to examine if patterns of collection Month and Location associated with sample ASV composition (Figure 3)
2. The Vegan function Adonis (multi factorial PERMANOVA) tested for Site and Month main effects for each distance estimates. Summaries are available for each measure in the file [adonis_data_allMetrics.csv](https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/text_tables/adonis/adonis_data_allMetrics.csv).
3. We also measured whether the significant main effects observed in the Adonis tests were due to dispersion, thus we also used the Vegan PERMDISP function to test if within-group distances to group centroid differ across groups. Site and Month effects were tested separately. Summaries are available in the file [betaDispersion_data_allMetrics.csv](https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/text_tables/permdisp/betaDispersion_data_allMetrics.csv).

# Machine Learning and Core Features
While the alpha and beta diversity tests provide some sense of how the community composition may vary in space and time, we wanted to better understand whether there was a smaller collection of sequence features that were more prevalant than others, and possibly, whether particular features were more associated with particular collection Sites or Months. There were two question we aimed to address:
  - First, what features are routinely present across samples. This was analyzed using QIIME's `core-features` function
  - Second, are there particular features that are more associated with a particular Site or Month? Because there are over 1,000 possible features, we restricted this analysis only to those samples present in at least 10% of all samples

## Core Features workflow 
Getting the core features in a dataset is simple with QIIME2's `core-features` function. However, becuase we were using a clustered dataset that was not rarefied, we first discarded any samples with less than 10,000 samples:
```
qiime feature-table filter-samples \
  --i-table /scratch/dro49/qiimetmp/mysotmp/Mangan.clust_p985_table_Filtd.qza \
  --p-min-frequency 10000 \
  --o-filtered-table Mangan.clust_p985_table_Filtd_min10k.qza
```

We then identify which sequence features are found among at least 10 percent of our samples. 
- First, identify core features among these remaining 196 samples
```
qiime feature-table core-features \
  --i-table Mangan.clust_p985_table_Filtd_min10k.qza \
  --p-min-fraction 0.1 \
  --p-steps 10 \
  --o-visualization Mangan.clust_p985_corefeatures.qzv
```

## Predicting sample metadata workflow 
- Second, export the visualization and capture those features present in at least 10% of samples
```
qiime tools export --input-path Mangan.clust_p985_corefeatures.qzv --output-path tmpdir_coreFeatures_min10k
cut -f 1 ./tmpdir_coreFeatures_min10k/core-features-0.100.tsv > featureIDs_core10.txt
```



- Lastly, filter those FeatureIDs from initial MachineLearn data to evaluate percent importance
  qiime feature-table filter-features \
  --i-table Mangan.clust_p985_table_Filtd_min10k.qza \
  --m-metadata-file featureIDs_core10.txt \
  --o-filtered-table Mangan.clust_p985_table_coreFeat10.qza
```

The QIIME2 [classify-samples](https://docs.qiime2.org/2020.8/tutorials/sample-classifier/) plugin offers a means to test a cross-validated supervised learning classifier - in this case, a Random Forest classifier. The goal is to take an input feature table and a set of metadata variables (in our case, the Site or Locations when guano was collected), train the classifier on some subset of these data, then test the classifier on the other subset of the data not used in training. The decision trees of the classifier used in the training set are then applied to the testing subset of data not used in training, and the accuracy of the classifier is assessed based upon how frequently the class variable prediction matches the truth. There are a few different flavors of cross validation available within the q2-sample-classifier within QIIME - we chose the nested cross validation (NCV) option to ensure that all feature data is used for both testing and training. This is possible because there are multiple iterations of testing and training - we used the default number of iterations (5) in our tests.

The individual features used in training the classifier are individually removed and the subsequent accuracy of the model lacking that information can be evaluated to the original model to determine each features relative importance to the predictive power of the classifier. This key feature enables us to determine which ASVs are more or less important at distinguishing among the class variables like Site or Month. Thus, while ordinations can help understand if the overall composition of samples differ by these variables, a supervised learning approach helps identify the individual ASVs that are driving the largest differences among these variables.

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
