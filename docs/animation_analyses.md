# Background
As explained in the [machine learning section](https://github.com/devonorourke/mysosoup/blob/master/docs/diversity_workflow.md#machine-learning) of the `diversity_workflow.md` document, the [machine_learn_analyses.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/machine_learn_analyses.R) script was used to generate a series of data summaries that highlighted the features (ASVs) that were most important to the model accuracy of the Random Forest classifier. We generated a series of static figures to illustrate changes in arthropod detections and abundances per Site and Month groupings (Figures [5](https://github.com/devonorourke/mysosoup/blob/master/figures/ml_pAbu_by_pOcr_scatterplot_byMonth_andSite.png) and [S4](https://github.com/devonorourke/mysosoup/blob/master/figures/ml_pAbu_by_pOcr_scatterplot_byMonth_andSite_DipteraOnly.png)), however these visualizations make it challenging to appreciate the magnitude of turnover among sequence variants in each group. Thus, we were motivated to create a few animations that make it easier to observe these trends.

# Animation-1: Month + Group, all Samples, all Orders in one plot
The ASVs highlighted in this first animation present the same dataset as in the static figure, but collapse the Month dimension in the animation (the Month in question is depicted in the title in the top left of the image). Each point represents a single ASV, with labels combining the associated Genus name (if available) as well as an ASV number that represents the relative rank of sequences that ASV contained across the entire dataset (i.e. ASV1 had the most sequences among all samples, ASV2 had the second most, etc.). We used solid points to represent the ASVs identified as important to the supervised learning model, while the square shapes represent ASVs that were not important to the model. While there were more than 8 arthropod Orders in the entire dataset, those selected represented the Orders with ASVs in the 75th percentile for model accuracy (importance).

![perSiteMonth_ASVchanges_MLonly](https://github.com/devonorourke/mysosoup/blob/master/figures/gifs/sitemonth_ASVs.gif)

There are several ways in which features are likely important for the classifier to make distinctions among groups:
1. Some features are highly abundant and frequently detected in all groups (e.g. `ASV-1|Eustala`), yet these proportions vary enough for differences to be meaningful for the decision tree to discriminate among Site and Month classes. For example, with `ASV-1|Eustala` the key distinction is likely the June to July transition between sites, where we see a _decrease_ in read abundance at the Egner site, but a slight increase at Hickory Bottoms.
2. Other groups tend to have Monthly changes in abundance or detections that follow similar patterns for each Site, but the magnitude of these differences are what help the model discriminate among Site and Month. For example, `ASV-4|Erioptera` tends to increase in abundance and detections from June into July and September at both Sites, however, the total number of detections and read abundances are generally larger for Hickory Bottoms than Egner.
3. Some features are simply unique to a given Site or Month. For instance, ASV's **12, 31, 33, and 47** are distinct to Egner, while ASVs **42, 56, 59, and 70** are private to Hickroy Bottoms. While some of these ASVs may be rare variants, many are not. Among those ASVs mentioned above:

| ASV | nSamples |
| --- | -------- |
| 12 | 28 |
| 31 | 1 |
| 33 | 4 |
| 37 | 8 |
| 42 | 11 |
| 56 | 1 |
| 59 | 4 |
| 70 | 4 |

# Animation-2: Month + Group, all Samples, each Order in unique facet
Because there are so many different ASVs in flux for a given Site or Month, we revised the animation above to subset each arthropod Order into it's own plot. Thus the same data used in the previous plot is used for the follow animation.

![perMonth_ASVchanges_perOrder](https://github.com/devonorourke/mysosoup/blob/master/figures/gifs/sitemonth_byOrder_ASVs.gif)

A few insights from this figure:
1. The Dipteran Order clearly have the most ASVs identified as important by the Random Forest classifier.
2. Nevertheless, other Orders indeed have taxa experience significant shifts in abundance and/or detections, though the diversity of these important ASVs are fewer per Order compared to Dipterans. It's also clear that the abundances or detections of these other Orders is relatively less than the Dipteran ASVs. Our previous work with mock communities suggest that the primers we're using are not biased towards Dipterans and amplify these Orders presented here (and others) equally well in vitro, thus these differences in abundance of sequences may be more due to biomass than molecular biases.


# Animation-3: Month + Group, all Samples, Diptera taxa only
To get a clearer perspective on the dynamics of Dipteran turnover, we generated the same plot as above but focused only on the Dipteran species. The labels on the left edge of the facet are assigned to taxa not identified by the classifier as important to the model, while the labels on the right edge of each facet were ASVs identified as relevant tot he classifier. As with previous plots, filled circles represent ASVs identified as relevant to the classifier, and open squares represent ASVs that did not contributing to the model importance.

![perMonth_ASVchanges_DipteranOnly](https://github.com/devonorourke/mysosoup/blob/master/figures/gifs/sitemonth_ASVs_DipteraOnly.gif)

As with the first animation, a similar set of patterns emerged:
1. Some ASV's have high abundances and detections in both Sites, but the directions of those changes differ by Month of collection. `ASV-3|Rhipidia` is more frequently detected and generates the most reads in Egner collections in June, while the same ASV is more abundant and detected in July for Hickory Bottom Samples.
2. Some ASVs like `ASV-2|unassigned` generate similar trends in detections but have distinct magnitudes with respect to read abundances. In July, for example, there were 12 samples collected at Egner with that ASV and 17 samples collected at Hickory Bottoms, yet Egner samples collectively generated 37,047 sequences compared to just 1,126 in Hickory Bottoms.  

Most ASVs did not contribute to model, and among the 1190 potential Dipteran ASVs, just 122 of these were used to train the classifier. This may partly be explained by the fact that only a subset of the samples were used for training, while other samples were excluded so that the testing would use independent datasets. However, it's also possible that some of the ASVs that were excluded simply did not have significant abundances or detections, and that is evident in the above figure. Among the more abundant or detected ASVs not identified as important by the classifier (e.g. `ASV-23|Helius`, `ASV-39|Nephrotoma`), we found that these ASVs tended to be in both the training set and the testing set, suggesting that their lack of importance may have to do with specific classifier decision parameters rather than just being left out from one of the two groups (training vs. testing).
