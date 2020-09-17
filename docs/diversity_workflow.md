# Background

Prior to beginning diversity analyses, we completed the following quality control steps:

- Raw sequence reads were trimmed with cutadapt and denoised with DADA2 into representative sequences resulting in the [Mangan.raw_linked_required.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.raw_linked_required.repSeqs.qza)) object, as described in the [sequence_processing.md](https://github.com/devonorourke/mysosoup/blob/master/docs/sequence_processing.md) workflow.    

- We identified bat-host sequences, as described in the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md) document. The host database design was described in the [host_database.md](https://github.com/devonorourke/mysosoup/blob/master/docs/host_database.md) document.  

- The non-bat representative sequences were investigated for contamination as described in the [contamination_investigations.md](https://github.com/devonorourke/mysosoup/blob/master/docs/contamination_investigations.md) document - we failed to detect extensive contamination either during DNA extraction or PCR amplification. See also the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script related to this contamination investigation.  

Given this information, we proceeded to filter the original/raw DADA2-processed feature table and sequence .qza artifacts as described in the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md#sequence-processing-of-non-bat-dna-prior-to-diversity-workflow) document:
- removed all batch and negative control samples, 
- clustered at 98.5%
- classified clustered representative sequences using a hybrid approach, retaining taxonomies prioritized by (1) exact alignments in VSEARCH, then (2) naive Bayes 
- filtered classified samples to retain only those samples assigned to the Arthropoda phylum, with at least family-level taxonomic names 
- rarefied samples using a depth of 10,000 sequences per sample. 

## Inputs for this workflow
- The rarefied table of clustered sequences, [Mangan.clust_p985_table_Rarefyd.qza](https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_table_Rarefyd.qza), is the file that is used for alpha and beta diversity estimates
- The non-rarefied table of clustered sequences, [Mangan.clust_p985_seqs_Filtd.qza](https://github.com/devonorourke/mysosoup/raw/master/data/qiime_qza/Mangan.clust_p985_seqs_Filtd.qza) is used for the `core-features` and `q2-classify-samples` analyses.


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
  --i-table Mangan.clust_p985_table_Filtd.qza \
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
