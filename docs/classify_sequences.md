# Filtering bat reads
We want to separate our amplicons into three groups:
1) Those amplicons that are likely to be associated to the bat host,
2) Those amplicons that are likely to be associated with the diet,
3) all amplicons that don't fit into 1 or 2.

To begin separating these data, we'll first try to identify which amplicons associate with host (bat) DNA. We created two separate databases - one primarily to identify bat host sequences, and another that contains millions more COI sequences from a variety of taxa including arthropods, chordates, fungi, and other eukaryotes. While there are overlaps in bat COI sequences between thos host-specific and the broad COI databases, we utilize the host-specific database as a first-pass to screen for likely host DNA sequences, then follow up with the larger COI database. The rationale is that each bat project may contain distinct hosts that may or may not be included in the broader database; rather than update the larger database for every project, we can screen out likely host sequences with this smaller database first.

## Filtering with host database
We begin by aligning all denoised ASVs against the host database containing a series of possible host reference sequences - see `host_database.md` for full details of how the samples were selected. These references included the focal species of _Myotis sodalis_ as well as a range of bats found across the Midwest and Northeast. In addition, we included other bat and bird species that had guano samples processed in our lab in previous projects; inclusion of these host sequences was done as a measure of contaminant check and removal. Samples that are at least 85% identical across at least 80% of the alignment subject are retained for analysis. Notably, this are _very_ liberal alignment parameters and will almost always include non-host sequences; however because the reference database contains only host-sequences, the classification output will suggest these are host samples. They aren't necessarily! However, we are casting as wide a net initially to retain any _possible_ ASVs that even resemble our host sequences, then follow up with more specific parameters (see below) to determine whether these sequences are truly host associated or not.

Representative ASVs from the `Mangan.samp-filt.repSeqs.qza` file were queried for alignment to host reference sequences:
> `$REF` refers to the file path to the relevant host files available at the _mysosoup_ GitHub repo

```
REF=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/host_dbs
REPSEQ=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.samp-filt.repSeqs.qza

qiime feature-classifier classify-consensus-vsearch \
--i-query "$REPSEQ" \
--i-reference-reads "$REF"/host_seqs.qza \
--i-reference-taxonomy "$REF"/host_taxonomy.qza \
--p-perc-identity 0.85 \
--p-query-cov 0.80 \
--p-maxaccepts 0 \
--p-strand both \
--p-threads 24 \
--o-classification host_tax_p85c80.qza
```

> Notably, we tested a range of percent identity parameters including 90%, 94%, and 98%. As percent identity is increased, the specificity of each ASV being classified is improved providng a better resolution of species; however fewer ASVs are retained. This information is useful for discerning the specific species of host DNA from each sample, and is used later once host-ASVs are identified and separated from non-host.

We export the `host_tax_p85c80.qza` file and search for any ASV which was retained _some_ level of classification (other than _Unassigned_); this list of ASVs is then used to filter the original fasta file of all DADA2-denoised representative sequences, and that "putative host" fasta file is then passed into NCBI BLAST to generate a list of candidate taxa with alignment percent identity and coverage values.

```
REPSEQDIR=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/mangan_link_req_repSeqs
TAXDIR=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/tax

cd $TAXDIR
qiime tools export --input-path host_tax_p85c80.qza --output-path host_tax_p85c80
grep -v 'Unassigned' ./host_tax_p85c80/taxonomy.tsv | cut -f 1 | tail -n +2 > potentialHits.txt
cd $REPSEQDIR
grep -f $TAXDIR/potentialHits.txt dna-sequences.fasta -A 1 | sed '/--/d' > $TAXDIR/potentialHits.fasta
```
There are 26 potential sequences identified in the resulting `Mangan_hosthits.fasta` file. These sequences were queried in NCBI BLAST; we retained the top 10 hits for each input in the `mangan_hittable.csv` file. As mentioned previously, it's possible that one or more of these ASVs are derived from non-bat COI sequences - possibly they are arthropod-derived sequences that share a relatively high degree of homology to the host-reference sequences. To distinguish the truly bat-associated sequences from non-bat sequences, we analyzed the taxonomies assigned to the best alignment parameters (coverage and percent identity ) using the `host_blastID.R` script. These analyses reveal that just 16 of 26 putative host ASVs are likely associated with bat hosts, and of those bat hosts, all three unique species had been previously detected in this area - see the full set of results available in the `Mangan.blastOut.csv` file. Notable, one and only one species is identified for each bat-associated ASV classified within NCBI; thus the 16 ASVs classified as some bat species are highly likely to be derived from a particular bat host. The remaining 10 ASVs represent a mix of likely arthropod-associated taxa, as well as likely eDNA from white-tail deer known to inhabit the area. The bat-associated ASVs are retained in the `host_hashIDs_ncbi.txt.gz` file.

With a list of likely bat associated ASVs (and other ASVs _not_ likely to be bat associated) we next explored two alternative classification approaches using a much larger database. The rationale is that our select (and sparse) host database may be insufficient and there may be other bat host sequences that could possibly be aligned. We utilized the same VSEARCH global alignment search program as done with the host database but using the broader COI references, and then explored a second classification approach (machine learning with a naive Bayes classifier) using the same broad reference database.


## Filtering with broad COI database using global alignment or naive Bayes classifier
While the initial approach using a host-associated database can help us filter out ASVs likely associated with a bat sample, the distinction between host and non-host is not equivalent to non-diet and diet components: there are likely several ASVs which are neither host nor diet-derived (for example, fungal ASVs and other eDNA which is amplified from passively collected guano). To better understand which ASVs are associated with host, with diet, or otherwise, we classified all denoised ASVs using two approaches: alignment-based with VSEARCH, and a naive-Bayes machine learning classifier. The database used in both of these methods was common; though this database represented a large collection of arthropod and non-arthropod samples. See `database_construction` for information on how the database used in classification was created.

### Alignment approach
All sequences were initially classified with VSEARCH, requiring a 97% identity and 94% coverage. These parameters were chosen because they reflect a reasonable tradeoff between specificity and breadth (making the identity too high will make reduce taxa ambiguity, but also discard samples less that fail to find a match at or above that identity threshold; making the identity too low will increase the number of samples with a match, but will also increase the ambiguity in the classification of each sample). Thus unlike the VSEARCH alignment using the host database, we are using more conservative alignment parameters when classifying all of our ASVs because we want to use these results to both distinguish bat host from non bat COI amplicons as well as use the resulting taxonomies to further filter and retain likely arthropod COI sequences (which we beleive to be associated with bat diets). Further filtering is explained later in this document.

```
READS=/mnt/lustre/macmaneslab/devon/guano/NAU/Mangan/qiime/reads/Mangan.samp-filt.repSeqs.qza
REFDIR=/mnt/lustre/macmaneslab/devon/guano/BOLDdb/bigDB

qiime feature-classifier classify-consensus-vsearch \
--i-query "$READS" \
--i-reference-reads "$REFDIR"/bigCOI.derep.seqs.qza \
--i-reference-taxonomy "$REFDIR"/bigCOI.derep.tax.qza \
--p-maxaccepts 0 \
--p-strand both \
--p-query-cov 0.94 \
--p-perc-identity 0.97 \
--p-threads 24 \
--o-classification Mangan_linked_tax_bigDB_p97c94.qza
```
> "$REFDIR" refers to reference file path to taxonomy mapping file and sequences
> "$READS" refers to path to `Mangan_noBats.repSeqs.qza` file


The resulting `.qza` artifact was exported and used in the `sequence_filtering.R` script for further processing.
```
qiime tools export --input-path Mangan_tax_p97c94.qza --output-path mangan_vs_linked_required
cd ./mangan_vs_linked_required/
mv taxonomy.tsv mangan_tax_vs.tsv
```

The most conservative parameters resulted in classifying just 5 ASVs from the dataset as some type of bat:
```
grep 'Chiroptera' mangan_tax_vs.tsv
```
These ASVs match the expected 3 bat species identified with the NCBI blast approach; furthermore all 5 ASVs in this match are identified in the previous classification approach using the alternative host database. However, we're missing 11 other expected bat ASVs - this is very likely a product of one of the two alignment parameters: the query coverage. Indeed, the 11 other ASVs classified as bat-associated from NCBI are "Unassigned" in this global alignment method. Given that the percent identities in these 11 ASVs from the NCBI approach were all above 97%, we suspect that the alignment length is the reason they aren't classified as such here; we can test this directly by running the same analysis again but lowering the query coverage parameter to the same value used in the earlier NCBI approach. We'll execute the same VSEARCH code as done previously, except we'll alter the query coverage parameter: `--p-query-cov 0.80`.

The resulting files illustrate that the query length is making a significant difference in the retention of the host fragments, but not altering the overall number of fragments classified overall:
  - for samples classified with a query coverage of 80% there are **1,494** *Unassigned* ASVs while **2,784** are *assigned*  
  - for samples classified with a query coverage of 94% there are **1,521** *Unassigned* ASVs while **2,757** are *assigned*
This amounts to just 27 different taxa being classified or not. However, while there were just 5 ASVs classified as a bat, the reduced query coverage classifier results in 14 ASVs classified as one of the three bat species we expect. Just two ASVs are missing from the NCBI approach (`9571a06c0df51336c4b5dae125329da3` and `076fa78f4dd69b67a085a99a0b560275`) while the other 14 (of 16) identified in the NCBI method were detected in with this distinct database. Interestingly, these two ASVs are _longer_ than the other fourteen bat ASVs: they're 223 bp instead of the expected 181 bp. Both of these sequences failed to have the 3'-end primer of the reverse complement reverse primer removed; apparently these few sequences failed to be properly removed from the Cutadapt script. These untrimmed sections drive the percent identity to be less than our 97% threshold (because 23 bases no longer align to the bat reference); indeed, by manually aligning those two remaining _Unclassified_ ASVs, they're clearly bat sequences like the other fourteen classified samples. This additional investigation supports the notion that there are 16 total ASVs expected to be classified as bats; however, to circumvent the issue of query coverage, we adopted a second classification approach that doesn't require specific alignment parameters: machine learning via a naive Bayes classifier.

### Machine learning approach

We also classified sequences using a Naive Bayes approach through the SciKit Learn package implemented in QIIME2. This classification strategy is fundamentally distinct from an alignment approach: the reference sequences are decomposed into discrete kmers and the resulting frequencies of these subreads are used to train the model. The representative sequences (our ASVs) are classified by associating the likelihood that those kmers are similar to the profiles of the reference sequences. Thus, while alignment confidence is a function of the actual length of coverage and percent identity between representative sequence and reference match, the confidence of the resulting machine learning classifier is a function of how well the kmer profiles of the representative sequences match those of the reference classifier.

Classification in this approach requires two steps: training the classifier (to build the models of taxonomy and kmer profiles), then applying the classifier to the representative sequences. To train the database using the combined arthropod, chordate, and host reference sequences:

```
qiime feature-classifier classify-sklearn \
--i-reads "$REPSEQ" --i-classifier "$CLSSFYR" \
--p-n-jobs 1 --p-reads-per-batch 2000 \
--o-classification mangan.sktax_linked_bigDB.qza
```

The trained classifier was then used to classify the representative sequences using default parameters:
> "$REPSEQ" refers to the full file path for the representative sequences `Mangan.samp-filt.repSeqs.qza`
> "$CLSSFYR" refers to the full file path for the trained classifier `nbClassifer_bigDB.qza`

The QIIME artifact was exported:
```
qiime tools export --input-path mangan.sktax_linked_bigDB.qza --output-path mangan_nb_linked_required
cd ./mangan_nb_linked_required/
mv taxonomy.tsv mangan_tax_nb.tsv
```

Unlike the default parameters for our VSEARCH approach, 16 and only 16 bat sequences are classified; in each case, these sequences are assigned to the species level with confidences exceeding 99.999%. A single non-bat mammal sequence is identified (the white-tail deer _Odocoileus virginianus_). These results mirror exactly what was discovered from NCBI. Collectively, these results suggest that there are 16 very likely bat-associated ASVs in the dataset. We next separate the bat ASVs from non-bat ASVs. Before moving into further filtering of non-bat sequences.

## Separating bat and non-bat data

Without our 16 bat ASVs identified we next isolate these ASVs to create a new ASV table and representative sequence file. These data are then used to understand which samples contain what bat species. We modify the `potentialHits.txt` file and use those 16 ASVs as a list to filter with. The file we create is called `bat_ASVs.txt` - note the addition of the `#OUTID` header:

```
#OTUID
076fa78f4dd69b67a085a99a0b560275
0f9b9d16fe3d116357d9fb68553da9e1
382ea10b1c328a71f139601c173f574b
414cdc433a260ea81a1df0130f7d054e
5b1b957e1f448c4404078ca20c5875c8
645527c3f16d27672ec42f6cb1252308
78659fdbc329b34d975ce914afbf78c6
9571a06c0df51336c4b5dae125329da3
9905172cd2a4905e4746823300fecfbf
b6d8cd13dd3c7ae8add173947e4453f0
ba3a67e00683bc23423185654463dbad
c1f82325df7f3777b7039ab8c6ebc15f
d2efef68ab6bafa1d39102b89bf88944
d522a80fed1adee71f631b85a3a16fed
de30f78335ddf188c9c3064447e5ecc9
ee6575cab7780fe25b51f8238ca74fd4
```

To create a bat-only ASV table and fasta files:
```
qiime feature-table filter-features \
--i-table Mangan.samp-filt.table.qza \
--m-metadata-file bat_ASVs.txt \
--o-filtered-table Mangan.batASVs.table.qza

qiime feature-table filter-seqs \
--i-data Mangan.samp-filt.repSeqs.qza \
--i-table Mangan.batASVs.table.qza \
--o-filtered-data Mangan.batASVs.repSeqs.qza
```

The same strategy is used to create a pair of files with the bat sequences **removed**:
```
qiime feature-table filter-features \
--p-exclude-ids \
--i-table Mangan.samp-filt.table.qza \
--m-metadata-file bat_ASVs.txt \
--o-filtered-table Mangan.nonbatASVs.table.qza

qiime feature-table filter-seqs \
--i-data Mangan.samp-filt.repSeqs.qza \
--i-table Mangan.nonbatASVs.table.qza \
--o-filtered-data Mangan.nonbatASVs.repSeqs.qza
```

# Next steps in analyses

1. The ` Mangan.batASVs*` files are used to determine which samples likely can be ascribed as a specific bat species - see the `host_seq-filtering.R` script. Executing this script produces a series of files:
  - `host_abundances-by-taxaHits.csv` represents the per-sample number of sequences of a bat-classified ASV (rows are samples, columns are bat species, elements of matrix are numbers of sequences); we discovered there are **200** instances with some bat ASV being present in a true sample (while 4 negative controls (of 7) contained some bat ASV also). This table further illustrates the instances in which multiple bat species may be assigned to a single sample; there are 18 such instances total. Notably, 5 of the 18 were from pooled samples (where multiple species isn't unexpected), and 2 others were from negative control samples. The remaining instances of multiple species assigned to a single pellet indicate are represented by one of two scenarios:
    A) there is a single dominant bat species with hundreds or thousands of reads for one bat species, and tens (or less) of sequences for another bat species
    B) there are less than 200 sequences in either sample
  The nature of sample collection (passive sampling where guano was piled in a mass beneath the roost) allows for such mutiple species COI detection to occur; nevertheless, the vast majority of samples were isolated from a single bat species, with 177 of 200 true samples being represented as _Myotis sodalis_ (though the range of sequence abundances therein ranged from just 3 to 27,774).
  - `sample_abundancesummaries_wBatHostData.csv` contains similar bat-species sequence counts, albeit with the addition of metadata (to indicate which samples came from which site and what kind of sample they were)


2. The `Mangan.nonbatASVs*` files are used in the R script `sequence_filtering.R` to determine which ASVs are considered diet components (arthropods) and which are likely non-diet (ex. fungal COI sequences, the deer, etc.). Notably, this file examines the extent of the ASVs associated with negative control samples and provides justification for why filtering these ASVs was likely unnecessary. A pair of text files containing lists of ASVs that are targets for our diversity estimates (of bat diet) are passed into the next part of the workflow: see `contamination_investigation.md` for details of how we used diversity tests to apply appropriate ASV filtering. The remaining filtered tables are then passed into `diversity_workflow.md` for complete analyses.