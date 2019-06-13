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
We began the contamination investigation by gathering some basic information about the negative control samples using the [sequence_filtering.R](https://github.com/devonorourke/mysosoup/blob/master/scripts/r_scripts/sequence_filtering.R) script. First, we calculated the number of unique ASVs in each sample and summed the read abundances of all ASVs per sample. As shown in the associated [data table](https://github.com/devonorourke/mysosoup/blob/master/data/text_tables/contam_evals/freqObs_and_seqDepth_perSample.csv) while all negative controls have total sequencing depths similar to many guano samples, the diversity of the negative controls (number of distinct ASVs) is vastly lower than most guano samples. In fact, there just two fo seven negative control sample that contains a high number of sequence variants (`ExtractionNTC1S1`and `ExtractionNTC8S64` contain 31 and 33 distinct ASVs, respectively), while all other negative controls contain 7 or fewer ASVs. This provides one piece of evidence against the hypothesis of pervasive contamination of DNA or PCR reagents, because we would not expect two negative control samples to have many more ASVs than other negative control samples. In addition, we determined that of all 87 ASVs present in negative control samples, there are never more than four negative control samples with a common ASV. In fact, just 4 of 87 ASVs are present in four samples, 5 of 87 ASVs are present in three samples, 15 of 87 are present in two samples, while 63 of 87 ASVs are detected in only one negative control sample. The broad uniqueness of sequence variants among the negative control samples is another strong indication that the detection of ASVs in negative controls is **not** a product of pervasive contamination, as we would have expected most ASVs in these samples to share similar sequence variants if that was the case.

We plotted the per-sample sequence abundance and number of distinct ASVs, but distinguished among those samples extracted by single tube extractions ("Isolate") versus those extracted by 96-well plate methods ("plate"). Among samples extracted by the plate method, we categorized wells that were surrounding a negative control well as a zone of contamination. Thus, any of the (up to eight) squares surrounding the central negative control well were marked as "TRUE" in the `ContamArea` legend below. In addition, we've altered the point type to clarify among guano samples from negative controls:

![contam_image](https://github.com/devonorourke/mysosoup/blob/master/figures/contam_eval_by_PlateType_and_ContamZone.png)





This workflow begins by creating a pair of filtered ASV tables:
1. True samples and control samples; remaining ASVs must have at least Family-rank taxa classification and

**edit below**
h the control samples and ASVs contained therein (as well as the true samples); a, then evaluates the appropriate sampling depth for each table. We then rarefy the datasets and move into the diversity assessments.
**edit above**

## Creating ASV-filtered tables and sequence files
First, generate the filtered tables using the text files output from the conclusion of the workflow described in the `classify_sequences.md` documentation:

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

Next, generate the remaining representative sequences:

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

These tables are then rarefied to an even sampling depth.


## Identifying appropriate sampling depths for rarefying

We'll first create a visualization to understand what sampling depth balances the need to preserve as many sequences as possible without discarding samples with less than some minimum sequence threshold:
```
qiime feature-table summarize --i-table Mangan.noNTCasvs-filt.table.qza --o-visualization Mangan.noNTCasvs-filt.sumry.qzv
qiime feature-table summarize --i-table Mangan.wNTCasvs-filt.table.qza --o-visualization Mangan.wNTCasvs-filt.sumry.qzv
```

We notice that the `Mangan.noNTCasvs-filt.table.qza` artifact has relatively fewer sequences and samples overall than the `Mangan.wNTCasvs-filt.table.qza` table (the table where all ASVs contained in any negative control sample were removed). It's likely that separate sampling depths are needed for rarefying these two tables. We want to maintain as many samples as possible without discarding diversity (ASVs). To explore the tradeoff between sampling depth and number of retained ASVs we use a QIIME script that generates an interactive alpha diversity rarefaction curve. By rarefying across a range of sampling depths we can determine how many ASVs are retained (depending on how many reads are subsampled). Because there are fewer reads in the `Mangan.noNTCasvs-filt.table.qza` table we reduced the difference between sampling depths for the visualization (note the different `--p-max-depth` parameters in the commands below):

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

The `Mangan.wNTCasvs-filt.alphaRareViz.qzv` visualization illustrates that a sampling depth of just 5,000 reads is likely sufficient to capture the majority of the diversity in every sample; however in the `Mangan.noNTCasvs-filt.table.qza` table we find that most samples retain the majority of their ASV diversity with as few as 1,000 reads. The lowered sampling depth in the later table ensures that we retain as many samples as possible while reducing the loss of diversity across the entire dataset; this makes comparing between the tables more direct because we have similar representation of samples per Site and Month.

## Rarefying filtered data tables
We rarefy each ASV table with the sampling depths determined above; note that in the case of the `Mangan.wNTCasvs-filt.table.qza` table we're creating two separate outputs: one with and one without any control samples (they are dropped from one of the tables). In the case of the `Mangan.noNTCasvs-filt.table.qza` table, all ASVs associated with negative control samples were removed thus these samples to not exist to be removed.

```
META=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/meta/qiime_meta.tsv

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
> The `$META` variable points to the `qiime_meta.tsv.gz` metadata file availabe in the metadata directroy of the GitHub repo


These tables serve as input for alpha and beta diversity estimates; they also are used in a machine learning process to determine which ASVs best predict a specified metadata factor (ex. Site or Date of collection). Each process is described below.

# Alpha diversity
Each filtered and rarefied `.qza` table is exported into an R environment where we apply the `Alpha-HillEstimates.R` script. In brief, this script calculates diversity estimates using three Hill Numbers (q=0, q=1, and q=2); higher q values represent increasing importance in relative abundances of reads per ASV in diversity estimates for a given sample. For q=0, the measure is equivalent to the number of observed ASVs. For q=1, the measure is equivalent to Shannon's index of diversity. For q=2, the measure is equivalent of Simpson's 1-D estimate of diversity. Importantly however, each q-value's metric can be directly compared - they share a common scale - as are estimators of the same (alpha) diversity metric.

We also explore whether the differences in alpha estimates vary by the factors of collection date (**Date**) and collection site (**Site**) using a a nonparametric Kruskal-Wallis test and a Dunn's test for pairwise significance testing, as well as a parametric multifactorial ANOVA. The overall impression is there is clear evidence for **Date** differences distinguishing among alpha diversity measures, though the degree that **Site** contributes to this change in diversity depends on the specific test:

Anova analyses:
  - These ANOVA suggests the observed differences in diversity are entirely attributed to the Date of collection, and have no Site differences
  - When ASVs present in NTCs remain there are Site and Date significant differences for q=1 and q=2, but not q=0. This suggests that there are some Site effects that are relevant to the ASVs we're removing with the filtering approach
  - These Site main effect is true whether we retain the negative control samples in the analysis or not.

Kruskal-Wallis analyses:
 - There are group differences for all q-levels for the 'noNTC' filtered data, with increasingly lower p-values as q-values are increased, so abundance information appears to matter in this context.
 - This abundance-related effect is REVERSED for the 'wNTC' data: including abundance information reduces group differences towards a non-significant threshold

Dunn pairwise comparisons:
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

We calculated distances directly in QIIME2 as well as using the Phyloseq package. Phylogenetic estimates for both packages required a rooted tree. Notably, the programs allow for trees to contain ASVs that are not present in the table, but all ASVs in the table must be present in the tree. Because the `Mangan.noNTCasvs.repSeqs.qza` table represents the source that the `Mangan.wNTCasvs.repSeqs.qza` table was filtered by, we create just a single:

```
REPSEQS=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.wNTCasvs.repSeqs.qza

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$REPSEQS" \
  --o-alignment aligned.NTCasvs.repSeqs.qza \
  --o-masked-alignment masked-aligned.wNTCasvs.repSeqs.qza \
  --o-tree unrooted-tree.wNTCasvs.qza \
  --o-rooted-tree rooted-tree.wNTCasvs.qza
```
> `$REPSEQS` represent the file path to the `Mangan.wNTCasvs.repSeqs.qza` file


The rooted tree file is exported as a Newick-formmated text file:
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
