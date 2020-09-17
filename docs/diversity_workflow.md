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
The QIIME2 [classify-samples](https://docs.qiime2.org/2020.8/tutorials/sample-classifier/) plugin offers a mechanism to perform cross-validated supervised learning - in our case, via a Random Forest classifier. The goal is to take an input feature table and a set of metadata variables (in our case, the Site or Month associated with a guano sample), and train a classifier on a subset of these data. The decision tree of the classifier model is then applied to the testing subset of data not used in training. While ordinations can help understand if the overall composition of samples differ by these variables, a supervised learning approach helps identify which sequence features are best associated with particular metadata variables.

There are a few different flavors of cross validation available within the q2-sample-classifier within QIIME - we chose the nested cross validation (NCV) option to ensure that all feature data is used for both testing and training. This is possible because there are multiple iterations of testing and training - we used the default number of iterations (5) in our tests. We are interested in two outputs of the sample classifier workflow: first, an evaluation of the overall accuracy of the classifier for a class of metadata; second, identifying which features are most important to building the classifier model.

We applied the same sample classifier workflow to three possible types of metadata:
- Site (Egner or Hickory)
- Month (June, July, or September) 
- Site + Month (Egner-June, Egner-July, Egner-September, Hickory-June, Hickory-July, Hickory-September)
> We did not further explore roost-specific characteristics as it was possible for multiple roosts to be occupied by the same group of bats.

To begin, the `core-features` visualization was exported to capture those features present in at least 10% of samples:
```
qiime tools export --input-path Mangan.clust_p985_corefeatures.qzv --output-path tmpdir_coreFeatures_min10k
cut -f 1 ./tmpdir_coreFeatures_min10k/core-features-0.100.tsv > featureIDs_core10.txt
```

We then filtered those FeatureIDs from our initial feature table to retain just those select features when predicting sample metadata using the q2-sample-classifier program:
```
  qiime feature-table filter-features \
  --i-table Mangan.clust_p985_table_Filtd_min10k.qza \
  --m-metadata-file featureIDs_core10.txt \
  --o-filtered-table Mangan.clust_p985_table_coreFeat10.qza
```

Lastly, we executed scripts in QIIME 2 using the default settings with one exception: we increased the number of estimators from 100 to 1000 to potentially increase our model accuracy by increasing the number of trees in each forest. Models for Month, Site, and SiteMonth, were performed similarly:
> substitute $METACOLUMN for one of three terms: 'SiteMonth', 'CollectionMonth', or 'Site'
> $METAFILE points to the [qiime_meta.tsv file](https://raw.githubusercontent.com/devonorourke/mysosoup/master/data/metadata/qiime_meta.tsv)

```
qiime sample-classifier classify-samples-ncv \
  --i-table Mangan.clust_p985_table_coreFeat10.qza \
  --m-metadata-file qiime_meta.tsv --m-metadata-column "$METACOLUMN" \
  --output-dir MLearn_ncv_"$METACOLUMN" \
  --p-parameter-tuning --p-estimator RandomForestClassifier --p-n-estimators 1000  

qiime tools export --input-path ./MLearn_ncv_"$METACOLUMN"/predictions.qza --output-path "$METACOLUMN"_Predictions
qiime tools export --input-path ./MLearn_ncv_"$METACOLUMN"/feature_importance.qza --output-path "$METACOLUMN"_featureImportance  
```

The six .tsv files exported to the three possible `"$METACOLUMN"_Predictions` and `"$METACOLUMN"_featureImportance` folters are are available at [the MachineLearn directory](https://github.com/devonorourke/mysosoup/tree/master/data/MachineLearn) of this repository. These data were used to generate the five panels in Figure 4:
- The three `*_importance.tsv` files were combined with metadata and read abundance data to create panels A and B in Figure 4 using the R script, [make_FeatureImportance_Heatmaps.R](https://raw.githubusercontent.com/devonorourke/mysosoup/master/scripts/r_scripts/make_FeatureImportance_Heatmaps.R).
- The three `*_predictions.tsv` files served as inputs to generating the heatmap plots shown in panels C-E of Figure 4 using the R script [make_Prediction_Heatmaps.R](https://raw.githubusercontent.com/devonorourke/mysosoup/master/scripts/r_scripts/make_Prediction_Heatmaps.R).
