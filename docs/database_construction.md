# Database curation
I've found that there are a handful of filtering steps that are useful when creating custom COI databases from BOLD records. To accomplish this, we used a few virtual environments for the work; `arrR` was used for data acquisition from BOLD, as well as filtering the resulting BOLD records,  `dev_qiime1` was used to apply further filtering parameters using some QIIME1 scripts and related Python scripts, and `qiime2-2018.11` was used to import the final filtered COI dataset. See the `sequence_filtering.md` document for QIIME2 information; the additional virtual environments were created as follows:

```
conda create -n arrR python=3.6
conda install -c r r-base
conda install -c bioconda seqkit
conda install -c bioconda csvtk
source activate arrR
```

However, we also required a few QIIME1 scripts as part of the dereplication process:
```
conda create -n dev_qiime1 python=2.7 qiime matplotlib=1.4.3 mock nose -c bioconda
```
> - I'd recommend testing proper installation, as I've encountered times where certain python packages need to be reinstalled
> - try launching the **dev_qiime1** environment, then type: `print_qiime_config.py -t`
> - if you get an error, uninstall and reinstall the offending python package (for example, I've had to do this for scipy and biom-format before)


## Obtaining data
We modified details illustrated in the vignette provided by the [bold R package](https://github.com/ropensci/bold) to download data from the [Barcode of Life Database](http://v4.boldsystems.org/). See the R scripts `bold_datapull.R` (Arthropods only), `bold_datapull_chordate.R` (Chordates only) and `bold_data_pull_nonArth_nonChord.R` (other animals, plus Fungi and Protist COI sequences) for full details on how the raw BOLD data was obtained and filtered. In brief, we selected COI records from BOLD, then data was filtered to require the `markercode` column matching "COI-5P". The output of this script produced `.csv` files that contain sequence information, taxonomic information, and associated metadata including the `sequenceID`, `processid`, `bin_uri`, `genbank_accession`, `country`, `institution_storing` records.

### Arthropod records
Because we're pulling nearly two million arthropod records from BOLD, we'll query the servers iteratively by generating a list of groups to pull from. I found that you can avoid any server complications by keeping the list size under about 2 million records, so I divvied up all Arthropods into the following groups:

1. all non-Insect arthropods
2. Dipterans
3. Coleopterans
4. Hymenopterans
5. Lepidopterans
6. all other Insects

See the R script `bold_datapull.R` for implementation on how data was queried; notably, I did not follow the vignette exactly. One particular sticking point was how the names are derived when using the `taxize` library - I found that it can erronously _miss_ some taxa because of different names being applied to a group. One such instance I identified occurred between the NCBI name of `Psocoptera` while BOLD uses the term `Psocodea`; the distinction is that NCBI things `Psocodea` is a superOrder instead of just an Order, and BOLD doesn't do the super/supra level distinctions. Because of this discrepency, you'd completely miss all `Psocodea` records if you followed the vignette exactly.

The retrieved records were split into two files: the first is a file containing sequence information with a particular format of taxonomic information, while the second contains metadata with taxonomy information in a different format. The first file taxonomic info is formatted in a way to produce the expected delimiters used with QIIME2 and Vsearch tools want (in case you want to dereplicate, for example), while the metadata file contained the various attributes for each Sequence ID, as well as the Phylum through Species taxa names in their original format (ie. not concatenated together). The latter allows us to quickly parse and tally the data as we impose certain filters.

## Filtering BOLD data
We filtered our raw BOLD sequence records by removing short sequences, sequences with gap characters, sequences with DNA characters that were not part of the standard IUPAC code, and sequences with extremely sparse taxonomic information.

### Filter1 - removing gaps and filtering for sequence length
The initial data set contains **3,175,054** records that passed our initial filtering requiring that a `markercode` column in the BOLD metadata be assgined with `COI-5P`. However, we identified a few hundred records in the initial BOLD data that contain gaps. We removed any gap character then selected only sequences longer than 100 bp to produce a single fasta file formatted for further processing in Vsearch. There were just **552** sequences that were discarded, so the vast  majority of the records (**3,174,502**) were retained:

```
grep -v 'sequenceID,taxon,nucleotides' boldCustom.allArth.seqNtaxa.csv | \
awk -F ',' '{OFS="\t"};{print $1";tax="$2, $3}' | \
sed 's/^/>/g' | tr '\t' '\n' | \
seqkit seq - --min-len 100 --remove-gaps -w 0 > boldCustom.allArth.filt1.fasta
```

### Filter2 - identifying any non-ATCG sequence characters
The BOLD sequences do not exlcusively contain non-ambiguous DNA characters. To see the distribution of which characters are present:
```
seqkit fx2tab -n -i -a boldCustom.allArth.filt1.fasta | \
csvtk -H -t grep -f 4 -r -i -p "[^ACGT]" | \
awk '{print $2}' | sort | uniq -c | sort -k1nr
```

About 90% of all records contain just ATCG characters, but the remaining contain either `N` or other ambiguous characters. Dig deeply into that list and you'll find that a single sequence contains a character that is not part of the standard DNA alphabet: `I` (it's the one with `ACGIT` characters). We need to remove that sequence entirely. To identify which sequence:
```
seqkit fx2tab -n -i -a boldCustom.allArth.filt1.fasta | \
csvtk -H -t grep -f 4 -r -i -p "I" | cut -f 1 > idsWithI.txt
```
> The result is record: `8982958;tax=k__Animalia;p__Arthropoda;c__Insecta;o__Diptera;f__Ephydridae;g__Hydrellia;s__Hydrellia			ACGIT`. This Order, Family, and Genus is well represented in the database, so dropping this one sequence is not going to affect our classification

We can discard that single sequence to create our fasta file filtered by length, gaps, and DNA alphabet:
```
seqkit grep --pattern-file idsWithI.txt --invert-match boldCustom.allArth.filt1.fasta > boldCustom.allArth.filt2.fasta
```

### Filter3 - removing sequences with taxonomic incompleteness
Prior to dereplicating data, we considered whether or not to impose an additional filtering parameter relating to taxonomic completeness - should we discard any sequences that contain no information at a particular taxonomic level. There are **3,175,053** total sequences in `boldCustom.allArth.filt2.fasta`, and while all of these sequences contain Phylum and Class information, **1,856** of these are missing Order information, and **122,709** were missing Family level information. The majority of sequences missing Order level information were associated with the _Arachnida_ Class (1,346), with _Collembola_ (281) and _Malacostraca_ (102) also missing several sequences.  Among the missing Family level data, the _Insecta_ class was missing the most sequences (62854), with the _Arachnida_ (29688), _Collembola_ (16749) and _Malacostraca_ (9622) again missing the most information among non-Insect arthropods. The Orders represented within the _Insecta_ Class with missing Family information were unequally distributed, with _Hempitera_ (35679) missing more than any other Order, though notably each Order has a unique number of total records.

We opted to create an filtered dataset to be used in dereplication in which all sequences that lacked Order or Family information were discarded. This ensured that dereplicating data with our lowest common ancestor (LCA) approach would not suffer from information loss in the situation where two identical sequences containing different taxonomic records - one without, say, Order level information, and another with complete information - was reduced solely because of different levels of information completeness (as opposed to distinct information at equivalent levels). Filtering follows a similar approach in identifying the sequence headers without Family info, then removing them from the fasta file:

```
seqkit fx2tab -n -i -a boldCustom.allArth.filt2.fasta | csvtk -H -t grep -f 1 -r -i -p "f__;" | cut -f 1 > idsNoFam.txt

seqkit -w 0 grep --pattern-file idsNoFam.txt --invert-match boldCustom.allArth.filt2.fasta > boldCustom.allArth.filt3.fasta
```

The resulting `boldCustom.allArth.filt3.fasta` now contains **3,051,814** Arthropod records with at least Family-level information, the proper DNA alphabet, and sequences longer than 100 bp (about 96% of what we started with).
> A `missingFamilyCounts.txt` documents the number of missing sequence records, grouped by taxonomic Orders. See the tidybug/data/databases

### Chordate records
The same principles were applied for Chordate BOLD records as with Arthropod records except that because the dataset was smaller we could pull the entire dataset from the BOLD API in a single batch. The order of operations:

1. Obtained all chordate records from BOLD API using the `bold_datapull_chordate.R` script (results in **215,031** records).

2. Filtered for length and removed gap characters (results in **215,028** records):
```
zgrep -v 'sequenceID,taxon,nucleotides' boldCustom.allChordate.seqNtaxa.csv.gz | \
awk -F ',' '{OFS="\t"};{print $1";tax="$2, $3}' | \
sed 's/^/>/g' | tr '\t' '\n' | \
seqkit seq - --min-len 100 --remove-gaps -w 0 | gzip > boldCustom.allChordate.filt1.fasta.gz
```

3. Identified any non-standard IUPAC DNA characters and removed such records (Inosine, "I" was discovered twice):
```
seqkit fx2tab -n -i -a boldCustom.allChordate.filt1.fasta.gz | \
csvtk -H -t grep -f 4 -r -i -p "[^ACGT]" | \
awk '{print $2}' | sort | uniq -c | sort -k1nr

seqkit fx2tab -n -i -a boldCustom.allChordate.filt1.fasta.gz | \
csvtk -H -t grep -f 4 -r -i -p "I" | cut -f 1 > idsWithI.txt

seqkit grep --pattern-file idsWithI.txt --invert-match boldCustom.allChordate.filt1.fasta.gz | gzip > boldCustom.allChordate.filt2.fasta.gz
```
> Of the two records, one contained many duplicate records of a shared Genus while the other contained no Family-level information and would be dropped by our downstream filtering requirements.

4. Removing sequences with taxonomic incompleteness (**188,779** remaining records):

```
seqkit fx2tab -n -i -a boldCustom.allChordate.filt2.fasta.gz | csvtk -H -t grep -f 1 -r -i -p "f__;" | cut -f 1 > idsNoFam.txt
seqkit -w 0 grep --pattern-file idsNoFam.txt --invert-match boldCustom.allChordate.filt2.fasta.gz | gzip > boldCustom.allChordate.filt3.fasta.gz
```

## Non Arthropod, Non Chordate records
The same principles were applied to BOLD records queried for non Arthropod/Chordate records, except we used a more liberal taxonomic rank requirement in filtering our data: instead of retaining only those with at least Family-rank information, we included all sequences assigned to at least Phylum. Because these data are used as filters to remove sequences we do not need to know more exclusive taxa Rank: provided they are a Mollusc, or a Fungi, for example, we know enough to  remove them from our dataset. Using a more liberal rank requirement allowed us to retain a broader set of representative sequences in these groups.

### For non-Animal records
```
zgrep -v 'sequenceID,taxon,nucleotides' boldCustom.allnonAnml.seqNtaxa.csv.gz | \
awk -F ',' '{OFS="\t"};{print $1";tax="$2, $3}' | \
sed 's/^/>/g' | tr '\t' '\n' | \
seqkit seq - --min-len 100 --remove-gaps -w 0 | gzip > boldCustom.allnonAnml.filt1.fasta.gz
```
> No additionanl steps were taken; no non-IUPAC characters were detected

There are just **6,436** of these records.

### For non-Arthropod/Chordate animal records

The initial **239,678** records are reduced to **239,651** records following length filtering, and then just a single additional record was dropped after filtering for IUPAC code.

```
zgrep -v 'sequenceID,taxon,nucleotides' boldCustom.allNonArthChordAnml.seqNtaxa.csv.gz | \
awk -F ',' '{OFS="\t"};{print $1";tax="$2, $3}' | \
sed 's/^/>/g' | tr '\t' '\n' | \
seqkit seq - --min-len 100 --remove-gaps -w 0 | gzip > boldCustom.allNonArthChordAnml.filt1.fasta.gz

seqkit fx2tab -n -i -a boldCustom.allNonArthChordAnml.filt1.fasta.gz | \
csvtk -H -t grep -f 4 -r -i -p "[^ACGT]" | \
awk '{print $2}' | sort | uniq -c | sort -k1nr > nonArthChordAnml.list

seqkit fx2tab -n -i -a boldCustom.allNonArthChordAnml.filt1.fasta.gz | \
csvtk -H -t grep -f 4 -r -i -p "I" | cut -f 1 > idsWithI.txt

seqkit grep --pattern-file idsWithI.txt --invert-match boldCustom.allNonArthChordAnml.filt1.fasta.gz -w 0 | gzip > boldCustom.allNonArthChordAnml.filt2.fasta.gz
```

This file was then subject to dereplication, as explained in the `Dereplicating with LCA` section below.

# Non-BOLD references from Porter Lab
Because there were so few COI records for non-Animals, we retained additional COI records from GenBank that were mined from Terri Porter's database - see [her workflow here](https://github.com/terrimporter/COI_NCBI_2018) and their CO1 mtDNA [dataset here](https://github.com/terrimporter/CO1Classifier). We used their v3.2 release [here](https://github.com/terrimporter/CO1Classifier/releases/tag/v3.2-ref). Notably, they required species-level information for their data, but the records were not dereplicated.

First, we gathered their data:

```
wget https://github.com/terrimporter/CO1Classifier/releases/download/v3.2-ref/CO1v3_2_training.tar.gz
tar xzf CO1v3_2_training.tar.gz
```

We then selected all non-Arthropod/Chordate records from the dataset to produce a fasta and taxonomy file used in subsequent dereplicating steps discussed below:
```
zcat mytrainseq.fasta.gz | paste - - | grep -v "Arthropoda" | grep -v "Chordata" | tr '\t' '\n' > tPorter.nonArth_nonChord.fasta

cat tPorter.nonArth_nonChord.fasta | paste - - | cut -f 1 -d ' ' | sed 's/>//' > left.tmp
cat tPorter.nonArth_nonChord.fasta | grep '^>' | cut -f 2- -d ' ' | cut -f 3- -d ';' > right.tmp
paste left.tmp right.tmp > tPorter.nonArth_nonChord.taxa

cat tPorter.nonArth_nonChord.fasta | grep -v "^>" > bottom.tmp
paste left.tmp bottom.tmp | sed 's/^/>/' | tr '\t' '\n' > tPorter.nonArth_nonChord.nolabel.fasta
```

This resulted in **177,427** additional COI records; though notably, there are (1) redundant sequences within the file, and (2) there may be redundant records with BOLD. These redundancies are removed following dereplication next, while also accounting for potentially shared taxa when common sequences contain distinct taxonomic information.

## Dereplicating with LCA
Databases can contain redundant sequences; dereplicating datasets is one solution to this problem, however, the default dereplication tools used in programs like VSEARCH will select the first record when multiple sequence matches exist. This can be problematic if the two records contain non-identical taxonomic information; this generally can create one of two problems. First, if the sequences contain non-identical records but equally complete levels of taxonomic information, there are two generally adopted strategies to picking a "winner" - majority, or consensus. A majority approach will take the most abundant of classifications, while a consensus approach invokes a least common ancestor (LCA) algorithm which retains only taxonomic information at the level where the matching sequences are equal. The second problem that can occur is if the two sequences contain unequal levels of information - a majority or consensus approach can again be invoked, but in this case you will always lose information with a consensus approach to whichever record contains the least amount of taxonomic information.
We adapted the [instructions](https://github.com/mikerobeson/Misc_Code/tree/master/SILVA_to_RDP) written by Mike Robeson for formatting a SILVA database, and incorporated the [consensus approach to reclassifying taxa](https://gist.github.com/walterst/9ddb926fece4b7c0e12c) script written by William (Tony) Walters. This resolved the problem of ties whereby identical sequences with distinct taxa are resolved to a common taxonomic level. This also required installing QIIME1 - we did so by creating a new virtual environment; see [here for instructions](http://qiime.org/install/install.html) on installing QIIME1. The actions are to first generate the appropriate file structures to perform the LCA algorithm on all data, then, with those "ties" now resolved, we can dereplicate the data knowing that any potential duplicate sequence selected will have the appropriate taxonomy assigned.

We created a single master file containing all COI records from all data sources by first creating the taxonomy mapping file and fasta files (required as input for the dereplication scripts) from Arthropod, Chordate, and other data sources separately (to speed up processing across multiple compute nodes). When combining all records we had to filter duplicate SequenceID records that arose from overlapping records obtained from our BOLD-derived data and the Porter Dataset that also pulled some of their records from BOLD (also detailed below) - there were just 120 instances of these among the 3.6 million records.

First, create the taxonomy mapping file; it's a two-column record of sequenceID and taxonomy string. The following code depicts this process for the Arthropod dataset only. The same process was repeated for Chordates, other Animal BOLD records, and the Porter dataset.

```
## Taxonomy mapping file:
cat boldCustom.allArth.filt3.fasta | grep '^>' | sed 's/^>//' | cut -d ';' -f 1 > tmp.left
cat boldCustom.allArth.filt3.fasta | grep '^>' | sed 's/^>//' | cut -d ';' -f2- | sed 's/tax=k__//' | sed 's/p__//' | sed 's/c__//' | sed 's/o__//' | sed 's/f__//' | sed 's/g__//' | sed 's/s__//' > tmp.right
paste -d '\t' tmp.left tmp.right | gzip > tmp_nolabels.taxa.gz
rm tmp.left tmp.right
```

Next, create a reduced fasta file without taxa info in the headers to save disk space when writing the subsequent files:
```
## reduced fasta file:
cat boldCustom.allArth.filt3.fasta | grep '^>' | cut -d ';' -f 1 > tmp.left
cat boldCustom.allArth.filt3.fasta | grep -v '^>' > tmp.seqs
paste -d '\t' tmp.left tmp.seqs | tr '\t' '\n' | gzip > tmp_nolabels.fasta.gz
rm tmp.left tmp.seqs
```

Once this process was completed for all subsets of files, we combined them into a single set of taxonomy map and fasta files:

```
cat tPorter.nonArth_nonChord.taxa.gz tmp_nolabels.taxa.gz boldCustom.allChordate_taxmap.txt.gz boldCustom.allNonArthChordAnml_taxmap.txt.gz boldCustom.allnonAnml_taxmap.txt.gz > allRecords_taxmap.txt.gz
```

Combine all fasta files:
```
cat tPorter.nonArth_nonChord.nolabel.fasta.gz tmp_nolabels.fasta.gz boldCustom.allChordate_nolabels.fasta.gz boldCustom.allNonArthChordAnml_nolabels.fasta.gz boldCustom.allnonAnml_nolabels.fasta.gz > allRecords_nolabels.fasta.gz
```

Results in **3,664,106** records. Unfortunately, a few of these are duplicates. To create deduplicated records for the fasta and taxonomy mapping file (these need to be decompressed for the  next process to function properly):
```
seqkit rmdup allRecords_nolabels.fasta.gz -w 0 > allRecords_nolabels_dedup.fasta
zcat allRecords_taxmap.txt.gz | sort -k1,1 | uniq > allRecords_taxmap_dedup.txt
```

This removed 120 duplicate records. Notably, duplicate sequences may remain - these deduplication here is with respect to the record's Sequence ID, not the sequence string itself. Before dereplicating, we need to identify duplicate sequence records and apply the LCA algorithm to reduce taxa information when applicable. This is a two step process that requires a pair of QIIME1 modified scripts; first, we use the `allRecords_nolabels_dedup.fasta` file as input for the [pick_otus.py](http://qiime.org/scripts/pick_otus.html) script. We switch the conda environments also:
```
conda activate dev_qiime1
pick_otus.py -i allRecords_nolabels_dedup.fasta -o pid100_otus --similarity 1.0 --threads 24
```

This generates a `allRecords_nolabels_dedup_otus.txt` text file within a newly created `pid100_otus` directory that is used in the next [create_consensus_taxonomy.py](https://gist.github.com/walterst/bd69a19e75748f79efeb) script to generate the mapping file that will apply the LCA algorithm.

We next apply three inputs (`allRecords_taxmap_dedup.txt`, `allRecords_nolabels_dedup`, and `allRecords_nolabels_dedup_otus.txt`) to generate a consensus mapping file output (`allRecords_taxmap_outmap.txt`) for our data:
```
python create_consensus_taxonomy.py allRecords_taxmap_dedup.txt allRecords_nolabels_dedup.fasta ./pid100_otus/allRecords_nolabels_dedup_otus.txt allRecords_outmap.txt
```

The `*outmap.txt`file is a 2-column file with the sequenceID and taxonomy string:
```
5333265	Animalia;Arthropoda;Insecta;Diptera;Cecidomyiidae;;
5333264	Animalia;Arthropoda;Insecta;Hemiptera;Aphididae;Aphis;Ambiguous_taxa
5333267	Animalia;Arthropoda;Insecta;Thysanoptera;Thripidae;;
```

The `allRecords_taxmap_dedup.txt` taxonomy mapping file (created from the (Record ID deduplicated) `allRecords_nolabels_dedup.fasta` file) contains **3,663,986** records - this includes potentially redundant sequences. The `*outmap.txt` file contains a similar number of records, but we have applied the LCA algorithm to these records. All that remains to do is dereplicate these records; because it doesn't matter which record is chosen among redundant sequences now (they all have similar taxa records) we can dereplicate the `allRecords_nolabels_dedup.fasta` directly:
```
vsearch \
--derep_fulllength allRecords_nolabels_dedup.fasta \
--output bigCOI.derep.fasta \
--relabel_keep --threads 4 --fasta_width 0 --notrunclabels
```

The dereplicated data contain **2,181,331** unique records (~ 60% of the original data); the Vsearch output provided the following brief summary:
```
2295494308 nt in 3663986 seqs, min 100, max 9861, avg 627
2181331 unique sequences, avg cluster 1.7, median 1, max 3547
```

## Formatting data for QIIME import
Two files are required for importing into QIIME2 to perform classification: (1) a taxonomy file, and (2) a fasta file. The taxonomy file uses the same 2 column format as the `*outmap.txt` file, except we're only going to retain the records present in the dereplicated dataset. The `bigCOI.derep.fasta` file contains the correct format for import.

We make a temporary list of all the sequenceID values from the headers of the `bigCOI.derep.fasta` file, then use that as a list to query the matches in the `*outmap.txt` file to generate the taxonomy file we want:
```
zgrep '^>' bigCOI.derep.fasta | sed 's/>//' > derep.seqid.tmp
fgrep -f derep.seqid.tmp allRecords_outmap.txt > bigCOI_taxmap.tmp
sed 's/Animalia;/k__Animalia;p__/' bigCOI_taxmap.tmp | sed 's/;/;c__/2' | sed 's/;/;o__/3' | sed 's/;/;f__/4' | sed 's/;/;g__/5' | sed 's/;/;s__/6' > bigCOI.derep_taxmap.txt
rm derep.seqid.tmp
```

I'm unclear exactly why there was a small discrepency between the number of records in the fasta and that in the taxonomy mapping file (six additional records in the tax map); these were identified and removed after manual inspection revealed their extra sub-species information was redundant:
```
## what's missing?
comm -23 <(sort taxmap.headers) <(sort fasta.headers) > extraTaxMap_data.txt

## produces:
2798360
2798361
2798362
2798372
2798374
2798376

## records are then dropped:
grep -v -f extraTaxMap_data.txt bigCOI_taxmap.txt > tmp && mv tmp bigCOI_taxmap.txt
```

> Those six duplicate records all had the taxa info: `Animalia;Arthropoda;Insecta;Trichoptera;Hydropsychidae;Hydropsyche;Hydropsyche nsp 2006031401`

One minor cleanup detail: we needed to convert the lower-case bases used in the Porter Dataset to upper case for an issue with classifying taxa:
```
seqkit seq bigCOI.derep.fasta --upper-case -w 0 > tmp && mv tmp bigCOI.derep.fasta
```
These `bigCOI.derep.fasta` and `bigCOI_taxmap.txt` files are imported into QIIME2:

```
## fasta file import
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path bigCOI.derep.fasta \
  --output-path bigCOI.derep.seqs.qza

## taxonomy file import
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path bigCOI_taxmap.txt \
  --output-path bigCOI.derep.tax.qza

## remove temporary file:
rm boldCOI.derep.qiime.tmp
```
