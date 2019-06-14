# Overview
It clear following the [classify_sequences.md](https://github.com/devonorourke/mysosoup/blob/master/docs/classify_sequences.md) workflow that some negative control samples contained sequence data that were shared among guano samples. However it wasn't clear whether these ASVs being detected in negative controls were indicative of pervasive contamination or simply isolated instances that did not impact the entire dataset. Though samples were processed using the same kit chemistry, they were extracted using two different methods: plate-based (via a 96-well kit) or isolated (via single tubes). All control samples were generated through the plate-based method.

There are a few different ways in which contamination may occur in our workflow, and as we highlight below, two of these are distinct to plate-based extractions:
- **PCR based**: If contamination was the result of some contaminant in the reagents prepping the amplicons (ex. contaminated primers, polymerases, BSA, etc.) we'd expect to see a high concentration of the same amplicons in every sample, because the same primers were used for all samples.
- **DNA based**: However, patterns of contamination via PCR may be distinct from contamination derived from DNA extraction methods, and importantly, there are multiple ways in which we may expect DNA-based contamination to be revealed. For instance,
  - **_reagent contamination_**: if reagents were contaminated during DNA extraction, we might expect all samples to contain similar abundances and taxonomic identities; unlike PCR however, we'd expect a variation in abundances of those high and low taxa, albeit the same taxa showing up repeatedly. However, this is unlikely because while we used the same kit chemistries, we did not use the same physical reagents (we first extracted by plate, then extracted by single tube). Thus, pervasive contamination would most likely be an indication of PCR issues if that was what the data suggested.
  - **_tip transfer contamination_**:if droplets of DNA were transferred with multichannel pipettes, then we'd expect samples in neighboring wells surrounding the negative control samples to share more ASVs in common than in non-neighboring wells
  - **_silicone mat contamination_**:if small chunks of guano were randomly dislodged from the silicone mat used in the plate-based extractions then we would expect certain negative controls to have higher abundances (those that had some guano fall into their empty well) than other control samples that didn't have any trace guano get dislodged. This issue of guano being dislodged can occur in the initial bead-beating step when the plate is shaken and the mat is removed to then begin the process of transferring supernatant into a new plate. In our experience even the most extensive centrifuging fails to completely pull down all the guano from the silicone mat.

With these patterns in mind, we used a series of measures to evaluate what to do about the negative controls. There were a few options:
1. Remove all ASVs detected in negative controls from the entire dataset.
2. Remove some number of reads for each ASV detected in a negative control from the dataset.
3. Keep all the ASVs found in both negative controls and guano samples.

## Evidence against pervasive contamination
We began the contamination investigation by gathering some basic information about the negative control samples using the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script. One ultimately useful data point was that we knew of a single example in which the technician performing the extractions watched a piece of guano fall into a control well. This was labeled as `blankS39`. Thus there were eight total negative control samples that generated sequence data, and one of those seven was _known_ to be contaminated.  

### Part 1 - most negative controls contain very low diversity of ASVs
First, we calculated the number of unique ASVs in each sample and summed the read abundances of all ASVs per sample. As shown in the associated [data table](https://github.com/devonorourke/mysosoup/blob/master/data/text_tables/contam_evals/freqObs_and_seqDepth_perSample.csv) while all negative controls have total sequencing depths similar to many guano samples, the diversity of the negative controls (number of distinct ASVs) is vastly lower than most guano samples. In fact, there just two of seven negative control sample that contains a high number of sequence variants (`ExtractionNTC1S1`and `ExtractionNTC8S64` contain 31 and 33 distinct ASVs, respectively), while all other negative controls contain 7 or fewer ASVs. Interestingly, the `blankS39` sample also contained a high number of ASVs (35), indicating that the two samples with higher numbers of unique ASVs may reflect instances in which small fragments of guano were dislodged from the silicone mat during the initial shaking step.
The stark differences in ASV counts (depsite similar sequence abundances) provides one piece of evidence against the hypothesis of pervasive contamination of DNA or PCR reagents, because we would not expect two negative control samples to have many more ASVs than other negative control samples. In addition, we determined that of all 87 ASVs present in negative control samples, there are never more than four negative control samples with a common ASV (including the `blankS39` sample). In fact, just 4 of 87 ASVs are present in four samples, 5 of 87 ASVs are present in three samples, 15 of 87 are present in two samples, while 63 of 87 ASVs are detected in only one negative control sample. The broad uniqueness of sequence variants among the negative control samples is another strong indication that the detection of ASVs in negative controls is **not** a product of pervasive contamination, as we would have expected most ASVs in these samples to share similar sequence variants if that was the case.


### Part 2 - zones of contamination don't associate with sequence diversity or read counts

We plotted the per-sample sequence abundance and number of distinct ASVs, but distinguished among those samples extracted by single tube extractions ("Isolate") versus those extracted by 96-well plate methods ("plate"). Among samples extracted by the plate method, we categorized wells that were surrounding a negative control well as a zone of contamination. Thus, any of the (up to eight) squares surrounding the central negative control well were marked as "TRUE" in the `ContamArea` legend below. In addition, we've altered the point type to clarify among guano samples from negative controls:

![contam_image](https://github.com/devonorourke/mysosoup/blob/master/figures/contam_eval_by_PlateType_and_ContamZone.png)

This plot provides some evidence suggesting against liquid handling contamination (i.e. errors by pipetting). There are no differences among guano samples within or away from negative control wells with respect to distinct sequence variants or sequencing counts; we'd expect the number of sequence variants among negative controls to be more  similar to those in their neighboring wells than those outside, but this is not observed. A better comparison would be to evaluate the actual sequence composition among these wells (which we do in a later section). Both of these lines of evidence suggest that there was not extensive pipetting errors for the samples extracted by the plate-based method.

Interestingly there is no differences in the distribution of per-sample numbers of reads (x-axis) or unique ASVs (y-axis) between extraction type (vertical facets, "Isolate" and "plate"). In fact, five of the eight negative controls have among the lowest numbers of sequence variants despite having moderate numbers of total reads. This is likely a result of two related items:
1. Pipette aerosols or liquids were unknowingly transferred into negative control samples through pipetting the various stages in which samples are incubated and plastic mats covering the plates are repeatedly attached and removed, or,
2. Trace guano fragments were dislodged from the silicone mat during the initial bead beating step

Trace amounts of contaminant solid guano or liquid extract could produce moderate sequence counts for a few reasons
  - Even small amounts of (contaminated) DNA in a negative control well can produce substantial number of amplicons, partly because the inhibitors present in guano are vastly diluted in the NTC well
  - DNA extracts were normalized prior to amplicifcation, and then PCR products were normalized again prior to pooling. Thus if PCR was even moderately successful for negative control samples and generated a relatively low concentration compared to most guano samples, these differences were erased prior to pooling.

In addition, the above plot highlights just three control samples with appreciable number of ASVs, and one of those three is the `blankS39` sample. It's therefore highly likely that those two additional negative control samples with similar diversity of sequence variants represent instances in which a fraction of a guano pellet was unintentionally transferred into a negative control well. This is not an indication of pervasive contamination, and the ASVs observed among those negative controls would not need to be removed from the entire dataset.


## Argument against dropping all ASVs detected in negative controls
At this point we were not concerned about pervasive contamination, however, we thought that perhaps we could just selectively remove the ASVs that are detected in the negative control samples. It turns out that's a terrible idea, because of the 68 ASVs detected in negative control samples _also present_ in guano samples, those ASVs accounts for over 45% of the entire library sequence abundance. In other words, most of the ASVs we're detecting in negative controls are also frequently detected in guano samples. Given that pervasive contamination is no longer likely, the idea of dropping the ASVs, or dropping some number of sequences from those ASVs would produce (needlessly) massive biases in our diet observations. What's much more likely is that these ASVs are detected in the negative control samples because a small fragment of guano indirectly fell into a well during DNA extraction, and because certain taxa (ASVs) are more likely to be consumed by a bat, they were more likely to be detected in these negative control samples.

You can see that clearly in the following plot showing the number of total sequences per ASV (x axis) and total number of samples and ASV is detected in (y axis). Points are colored as follows:
- red dots represent ASVs detected only in the eight negative controls (they generate virtually no reads, and are detected in fewer than 4 samples each)
- blue dots represent ASVs detected only in guano samples
- black dots are common to both guano and negative control samples

![contam_image2](https://github.com/devonorourke/mysosoup/blob/master/figures/contam_eval_ASV_and_SeqCounts_uniqnessColored_per_Guano_or_Control.png)

The common ASVs are routinely those that are present in high abundances and in many samples.

Because our downstream filtering for diet-taxa will require at least Family-level analyses, we concluded the sequence filtering R script by generating a pair of files:
- [taxfiltd_ASVs_NTCdrops.txt](https://github.com/devonorourke/mysosoup/blob/master/data/taxonomy/taxfiltd_ASVs_NTCdrops.txt): the ASVs that meet our taxonomy-filtering criteria, but any ASV detected in an NTC is removed
- [taxfiltd_ASVs_NTCincluded.txt](https://github.com/devonorourke/mysosoup/blob/master/data/taxonomy/taxfiltd_ASVs_NTCincluded.txt): the ASVs that meet our taxonomy-filtering criteria, but any ASV detected in an NTC is removed

We used those two files in the following diversity measures to provide an additional means with which we could evaluate the appropriate way to filter these ASVs.

## part 3 - diversity measures fail to detect trends in taxonomic composition
As mentioned earlier, we can use diversity estimates (both alpha and beta) to compare the negative control samples to guano samples. For instance, by calculating dissimilarity (distances) among all samples and ordinating, we can evaluate whether or not the community composition of samples surrounding the wells of the negative control samples are more similar to the negative controls (an indication of local contamination in which we might just drop certain samples), or not. We investigated estimates of alpha and beta diversity below using a series of metrics to better understand whether these ASVs shared among negative control samples and guano samples should be dropped or not.

### Creating ASV-filtered tables and sequence files
First, we filtered the raw QIIME-formatted sequence and table artifacts using the (`taxfiltd_ASVs_NTC*.txt`) text files output from the workflow above:

First the tables:
```
qiime feature-table filter-features \
--i-table Mangan.nonbatASVs.table.qza \
--m-metadata-file taxfiltd_ASVs_NTCdrops.txt \
--o-filtered-table Mangan.noNTCasvs-filt.table.qza

qiime feature-table filter-features \
--i-table Mangan.nonbatASVs.table.qza \
--m-metadata-file taxfiltd_ASVs_NTCincluded.txt \
--o-filtered-table Mangan.wNTCasvs-filt.table.qza
```

Then the classify_sequences
```
qiime feature-table filter-seqs \
--i-data Mangan.nonbatASVs.repSeqs.qza \
--i-table Mangan.noNTCasvs-filt.table.qza \
--o-filtered-data Mangan.noNTCasvs.repSeqs.qza

qiime feature-table filter-seqs \
--i-data Mangan.nonbatASVs.repSeqs.qza \
--i-table Mangan.wNTCasvs-filt.table.qza \
--o-filtered-data Mangan.wNTCasvs.repSeqs.qza
```

The resulting artifacts now have been filtered according to our specific taxonomic requirements (must be associated with the Arthropoda Phylum, must have at least Family-name information). The `wNTCasvs.*.qza` artifacts

Prior to calculating diversity estimates we first determiend appropriate sampling depths for each dataset.


### Identifying appropriate sampling depths for rarefying

We used a QIIME 2 program to generate a species accumulation curve (in this case, and ASV accumulation curve) to determine what sampling depth balances the need to preserve as many sequences as possible without discarding samples with less than some minimum sequence threshold. This program allows you to define the intervals of the subsampling depth, so we first generated a pair of summary visualizations in QIIME 2 that have an interactive component where you can quickly determine how many samples remain at a defined read depth. The summary visualizations were created by executing this code:
```
qiime feature-table summarize --i-table Mangan.noNTCasvs-filt.table.qza --o-visualization Mangan.noNTCasvs-filt.sumry.qzv
qiime feature-table summarize --i-table Mangan.wNTCasvs-filt.table.qza --o-visualization Mangan.wNTCasvs-filt.sumry.qzv
```

It's evident that the [Mangan.noNTCasvs-filt.table.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.noNTCasvs-filt.table.qza) associated with the [Mangan.noNTCasvs-filt.table.qzv](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qzv/Mangan.noNTCasvs-filt.sumry.qzv) artifact has relatively fewer sequences and samples overall than the [Mangan.noNTCasvs-filt.table.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table.qza) table (compare to the [Mangan.wNTCasvs-filt.table.qzv](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qzv/Mangan.wNTCasvs-filt.sumry.qzv) file, i.e. the table where all ASVs contained in any negative control sample were removed). That's expected, given that we discarded many ASVs that we know are present in many samples.

To explore the tradeoff between sampling depth and number of retained ASVs we use a QIIME script that generates an interactive alpha diversity rarefaction curve. By rarefying across a range of sampling depths we can determine how many ASVs are retained (depending on how many reads are subsampled). Because there are fewer reads in the `Mangan.noNTCasvs-filt.table.qza` table we reduced the difference between sampling depths for the visualization (note the different `--p-max-depth` parameters in the commands below):

```
qiime diversity alpha-rarefaction \
--p-metrics observed_otus \
--p-min-depth 100 --p-max-depth 5000 \
--i-table Mangan.noNTCasvs-filt.table.qza \
--o-visualization Mangan.noNTCasvs-filt.alphaRareViz.qzv

qiime diversity alpha-rarefaction \
--p-metrics observed_otus \
--p-min-depth 100 --p-max-depth 10000 \
--i-table Mangan.wNTCasvs-filt.table.qza \
--o-visualization Mangan.wNTCasvs-filt.alphaRareViz.qzv
```

The [Mangan.wNTCasvs-filt.alphaRareViz.qzv](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qzv/Mangan.wNTCasvs-filt.alphaRareViz.qzv) visualization illustrates that a sampling depth of just 5,000 reads is likely sufficient to capture the majority of the diversity in every sample associated with the `Mangan.wNTC` dataset; however the [Mangan.noNTCasvs-filt.table.qzv](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qzv/Mangan.noNTCasvs-filt.alphaRareViz.qzv) visualization suggests that the majority of ASV diversity is captured with as few as 1,000 reads. The lowered sampling depth in the later table ensures that we retain a similar number of samples between the two datasets while minimizing the loss of diversity. Thus we ensure a similar representation of numbers of samples per Site and Month when conducting the diversity estimates.

## Rarefying filtered data tables
We rarefy each ASV table with the sampling depths determined above; note that in the case of the `Mangan.wNTCasvs-filt.table.qza` table we're creating two separate outputs: one with and one without any control samples (they are dropped from one of the tables). In the case of the `Mangan.noNTCasvs-filt.table.qza` table, all ASVs associated with negative control samples were removed thus these samples to not exist to be removed.

> `$META` represents [qiime_meta.tsv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/qiime_meta.tsv), a QIIME-formatted metadata file similar to the original [mangan_metadata.csv](https://github.com/devonorourke/mysosoup/blob/master/data/metadata/mangan_metadata.csv) file.

```
qiime feature-table rarefy \
--i-table Mangan.noNTCasvs-filt.table.qza \
--p-sampling-depth 1027 \
--o-rarefied-table  Mangan.noNTCasvs-filt.rarefied-table.qza

qiime feature-table rarefy \
--i-table Mangan.wNTCasvs-filt.table.qza \
--p-sampling-depth 5200 \
--o-rarefied-table  Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza

qiime feature-table filter-samples \
--i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
--m-metadata-file $META \
--p-where "SampleType='control'" \
--p-exclude-ids \
--o-filtered-table Mangan.wNTCasvs-filt.rarefied-table_noNegSamps.qza
```

These tables serve as input for alpha and beta diversity estimates.

# Alpha diversity
Each filtered and rarefied `.qza` table is exported into an R environment where we apply the [contam_workflow_alphadiv.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/contam_workflow_alphadiv.R) script. In brief, this script calculates diversity estimates using three Hill Numbers (q=0, q=1, and q=2); higher q values represent increasing importance in relative abundances of reads per ASV in diversity estimates for a given sample. For q=0, the measure is equivalent to the number of observed ASVs. For q=1, the measure is equivalent to Shannon's index of diversity. For q=2, the measure is equivalent of Simpson's 1-D estimate of diversity. Importantly however, each q-value's metric can be directly compared - they share a common scale - as are estimators of the same (alpha) diversity metric.

![alphaHill_contam3](https://github.com/devonorourke/mysosoup/blob/master/figures/taxfilt_Alpha_Hillvals_wNTCasvs.png)

The plot above shows diversity estimates (y axis) for each of the three Hill Numbers (vertical facets)) for each sample, binned according to the roost it was obtained (x axis), and colored according to the month it was collected. The control samples generally lower diversity estimates compared to most guano samples for each of the Hill Numbers. The exception to this concerns two control samples for `q=1` and `q=2` Hill Numbers. These two samples are `blankS39` (the negative control we _know_ was contaminated with a piece of guano), and `ExtractionNTC8S64`, one of the two other negative control samples that had higher numbers of ASVs than the remaining controls. The higher diversity estimates with higher Hill values is a strong indication that this control sample is likely contaminated by an unintended guano fragment that was dislodged by a silicone mat during DNA extraction. Likewise, the relatively low diversity estimates for other negative control samples (compared to guano samples) suggests that pervasive contamination is not likely.

---

We also explored whether the differences in alpha estimates vary by the factors of collection date (**Date**) and collection site (**Site**) using a a nonparametric Kruskal-Wallis test and a Dunn's test for pairwise significance testing, as well as a parametric multifactorial ANOVA. The overall impression is there is clear evidence for **Date** differences distinguishing among alpha diversity measures, though the degree that **Site** contributes to this change in diversity depends on the specific test:

**Anova analyses**: Summaries available [here](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/anovas/contam_checks)
  - These ANOVA suggests the observed differences in diversity are entirely attributed to the Date of collection, and have no Site differences
  - When ASVs present in NTCs remain there are Site and Date significant differences for q=1 and q=2, but not q=0. This suggests that there are some Site effects that are relevant to the ASVs we're removing with the filtering approach
  - These Site main effect is true whether we retain the negative control samples in the analysis or not.

**Kruskal-Wallis analyses:** Summaries available [here](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/kruskal/contam_checks)
 - There are group differences for all q-levels for the 'noNTC' filtered data, with increasingly lower p-values as q-values are increased, so abundance information appears to matter in this context.
 - This abundance-related effect is REVERSED for the 'wNTC' data: including abundance information reduces group differences towards a non-significant threshold

**Dunn pairwise comparisons:** Summaries available [here](https://github.com/devonorourke/mysosoup/tree/master/data/text_tables/dunn/contam_check)
  - The majority of significant pairwise comparisons exist between Site and Dates that are not overlapping on a calendar: June and September. There are differences in diversity estimates across a range of Hill Numbers, and these differences tend to occur mainly for both sites between the June and September months (that is, there are significant differences between those Months within and between Sites). The number of differences are reduced among the `wNTC` dataset, suggesting that including the ASVs present in negative controls may drown out the signal. This doesn't necessarily indicate these are especially contaminated features - simply that they are likely present in so many samples and at such a high abundance that their inclusion is swamping the ability to detect meaningful differences among Site and Date groups.

# Beta diversity
We assessed beta diversity using four different measures: two non-phylogenetic and three phylogenetic.
  - Non-phylogenetic included:
    - Dice-Sorensen (**DS**): an occurrence-based metric that is equivalent in calculation to Bray Curtis without including abundance information
    - Bray-Curtis (**BC**): equivalent in calculation to **DS** except values are weighted by abundance information
  - Phylogenetic measures that incorporate branch lengths into diversity estimates relied on Unifrac, using:
    - Weighted unifrac (**WU**) included abundance information in diversity estimate
    - Unweighted unifrac (**UU**) occurrence data only (no abundance information)
    - Generalized unifrac (**GU**) includes abundance information but avoid weighting too much to rare or abundant lineages

We calculated distances directly in QIIME2 as well as using the Phyloseq package. Phylogenetic estimates for both packages required a rooted tree. Notably, the programs allow for trees to contain ASVs that are not present in the table, but all ASVs in the table must be present in the tree. Because the [Mangan.noNTCasvs.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/asvTables/Mangan.wNTCasvs-filt.table.qza) table represents the source that the `Mangan.wNTCasvs.repSeqs.qza` table was filtered by, we create just a single tree:

> `$REPSEQS` represents the [Mangan.wNTCasvs.repSeqs.qza](https://github.com/devonorourke/mysosoup/blob/master/data/qiime_qza/seqs/Mangan.wNTCasvs.repSeqs.qza) file containing all taxonomy-filtered non-bat sequence variants

```
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$REPSEQS" \
  --o-alignment aligned.NTCasvs.repSeqs.qza \
  --o-masked-alignment masked-aligned.wNTCasvs.repSeqs.qza \
  --o-tree unrooted-tree.wNTCasvs.qza \
  --o-rooted-tree rooted-tree.wNTCasvs.qza
```

The rooted QIIME-formatted rooted tree file is exported as [rooted-tree.wNTCasvs.nwk](https://github.com/devonorourke/mysosoup/blob/master/data/trees/rooted-tree.wNTCasvs.nwk), a Newick-formmated text file:
```
qiime tools export --input-path rooted-tree.wNTCasvs.qza --output-path rooted-tree.wNTCasvs.nwk
```

# Generating data for Adonis
We used distances estimated in QIIME2 for our Adonis evaluations. These included both phylogenetic and non-phylogenetic estimates.

## Non phylogenetic estimates of beta diversity in QIIME2
DS and BC distance measures were estimated for both `wNTC` and `noNTC`-filtered datasets in QIIME2; unfortunately there is no Morisita-Horn method available in QIIME2 (but there is in the Vegan-dependent Phyloseq R package):
```
qiime diversity beta \
--i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
--p-metric dice \
--o-distance-matrix Mangan.wNTCasvs.ds.beta.qza

qiime diversity beta \
--i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
--p-metric braycurtis \
--o-distance-matrix Mangan.wNTCasvs.bc.beta.qza

qiime diversity beta \
--i-table Mangan.noNTCasvs-filt.rarefied-table.qza \
--p-metric dice \
--o-distance-matrix Mangan.noNTCasvs.ds.beta.qza

qiime diversity beta \
--i-table Mangan.noNTCasvs-filt.rarefied-table.qza \
--p-metric braycurtis \
--o-distance-matrix Mangan.noNTCasvs.bc.beta.qza
```

## Phylogenetic estimates of beta diversityin QIIME2
Weighted and unweighted unifrac assesments used in these analyses:

```
TREE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/diversity/tree/rooted-tree.wNTCasvs.qza

qiime diversity beta-phylogenetic \
--i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
--i-phylogeny "$TREE" \
--p-metric unweighted_unifrac \
--o-distance-matrix Mangan.wNTCasvs.uu.beta.qza

qiime diversity beta-phylogenetic \
--i-table Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza \
--i-phylogeny "$TREE" \
--p-metric weighted_unifrac \
--o-distance-matrix Mangan.wNTCasvs.wu.beta.qza

qiime diversity beta-phylogenetic \
--i-table Mangan.noNTCasvs-filt.rarefied-table.qza \
--i-phylogeny "$TREE" \
--p-metric unweighted_unifrac \
--o-distance-matrix Mangan.noNTCasvs.uu.beta.qza

qiime diversity beta-phylogenetic \
--i-table Mangan.noNTCasvs-filt.rarefied-table.qza \
--i-phylogeny "$TREE" \
--p-metric weighted_unifrac \
--o-distance-matrix Mangan.noNTCasvs.wu.beta.qza
```
> The `$TREE` file represents the full path to the rooted tree generated previously (`rooted-tree.wNTCasvs.qza`)

All non-phylogenetic and phylogenetic distance estimates (as `.qza` artifacts) were used in `adonis_estimates.R` script.

# NMDS plots with phyloseq
We used the Phyloseq library in an R environment to calculate distances from the same `Mangan.wNTCasvs-filt.rarefied-table_wNegSamps.qza` and `Mangan.noNTCasvs-filt.rarefied-table.qza` ASV artifact files as input. See the `nmds_plots.R` script for full details.


# Conclusions of contamination findings
1. some controls have very few numbers of reads
2. some controls have average read abundance, but very few ASVs
3. some controls have average reads and average ASVs
4. ordination of all samples shows association of Month regardless of whether or not contaminated ASVs are included or not
5. ordinating by plotting well location shows no association to neighboring samples - it's not pervasive dropps of liquid
6. ordinatinog plots with extraction type (isolate or plate) shows no association with control samples - so it's not pervasive plate contamination

conclusions: likely two things happening
A. aerosols with low extract inhibition; note that samples were normalized
B. chunks of poop during plate loading, mat removal and transfers




## Machine learning

```
------
#!/bin/bash

#SBATCH -D /mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/classify-samples
#SBATCH --job-name="lrnBatch"

module purge
module load anaconda/colsa
source activate qiime2-2019.1


READPATH=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads
METAFILE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/meta/qiime_meta.tsv

qiime sample-classifier classify-samples \
  --i-table "$READPATH"/Mangan_noBats_famOnly_min9000seq.repSeqs.qza \
  --m-metadata-file "$METAFILE" \
  --m-metadata-column BatchType \
  --p-optimize-feature-selection \
  --p-parameter-tuning \
  --p-estimator RandomForestClassifier \
  --p-n-estimators 1000 \
  --output-dir learn-BatchType
```



Pulling out the sample names from the .qza file:
```
qiime tools export --input-path Mangan_noBats_famOnly_min9000seq_noControls.repSeqs.qza --output-path Mangan_noBats_famOnly_min9000seqs_noControls

biom convert -i feature-table.biom --to-tsv -o Mangan_noBats_famOnly_min9000seqs_noControls.tsv

head -2 Mangan_noBats_famOnly_min9000seqs_noControls.tsv | tail | tr '\t' '\n' | tail -n +3 > Mangan_noBats_famOnly_min9000seqs_noControls_sampleNames.txt
```



```




Estimate distances with ASVtable and rooted tree using weighted and unweighted unifrac:

TABLE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan_noBats_ASVtable_rarefied.qza
TREE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/diversity/tree/rooted-tree.qza

```
qiime diversity beta-phylogenetic \
--i-table "$TABLE" \
--i-phylogeny "$TREE" \
--p-metric weighted_unifrac \
--o-distance-matrix mangan_nobats_wuni_dist.qza

qiime diversity beta-phylogenetic \
--i-table "$TABLE" \
--i-phylogeny "$TREE" \
--p-metric unweighted_unifrac \
--o-distance-matrix mangan_nobats_uuni_dist.qza
```





/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/diversity/beta

unweighted_unifrac

QUESTIon:
I need to filter out the couple of mule deer reads before really doing the tree alignment and diversity estimates... those have to be removed from the rep-seqs and rep-table file...

Those are just 2 of over 5000 ASVs.
How should I go about removing any others?
Is it worth dropping the ASVs that are no longer classified by either alignment or Bayesian?
Could we use a tree built from ALL ASVs to identify the cluster that are clearly from arthropods against those that are likely NOT?





TABLE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan_noBats_ASVtable_FamOnly_rarefied.qza
TREE=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/diversity/tree/ASVfamOnlyTree/FamOnly.rooted-tree.qza

qiime diversity beta-phylogenetic \
--i-table "$TABLE" \
--i-phylogeny "$TREE" \
--p-metric weighted_unifrac \
--o-distance-matrix mangan_nobats_wuni_dist.qza

qiime diversity beta-phylogenetic \
--i-table "$TABLE" \
--i-phylogeny "$TREE" \
--p-metric unweighted_unifrac \
--o-distance-matrix mangan_nobats_uuni_dist.qza
