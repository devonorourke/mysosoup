# Host COI data
We included fasta references for bat and bird COI sequences from a variety of sources. These data included:
- Bat data from
  - Clare, Lim, Fenton, and Hebert 2011 PLoS supplementary list: CACA, CAPE, CASU, DERO, EPBR, EPFU, EUMA, GLSO, LABL, MOMO, MYNI, MYRI, NOAL, NOLE, PHDI, PHHA, PHST, PTGY, PTPA, PTPE, RHPU, RHNA, RHTU, STLI, TOSA, TRCI
  - from PopSet 726974368: MYSO, MYSE, and MYLE
  - from PopSet 301344216: MYAU, MYLU, PESU, EPFU, NYHU, LABO, COTO, LACI, LANO
  - from Chambers et al. (2015) - Species from Feces: EUGL, EUUN, RHBI, STPA
  - from Korstian, Williams, and Hale: LAEG
  - from BOLD (direct submission): MOSI

Additional host COI reference sequences for a few bird species were added because samples processed in the same lab as the Mangan study:
  - From Kerr, 2007: ABTO, ATFL, LUWA
  - From Kerr (2009): CORF
  - From Kerr, 2011: MAWA, JAWE, SOSP, WCSP, WIWA, YEWA
  - From Hebert 2004: BEVI,
  - Schindel 2011 (bioproject PRJNA81595): BEWR, YBCH
  - Direct BOLD submission: BHCO, BHGR, BLPH, CLSW, COYE, HAAM, LEGO, MELT, RBLE, VERD
  - From Saitoh et al. (2015: JABW  

Species codes used above refer to the following species:

| code | taxa name |
| --- | --- |
|ABTO|Melozone aberti|
|ATFL|Myiarchus cinerascens|
|BEVI|Vireo bellii|
|BEWR|Thryomanes bewickii|
|BHCO|Molothrus ater|
|BHGR|Pheucticus melanocephalus|
|BLPH|Sayornis nigricans|
|CACA|Carollia castanea|
|CAPE|Carollia perspicillata|
|CASU|Carollia subrufa|
|CLSW|Petrochelidon pyrrhonota|
|CORF|Carpodacus erythrinus|
|COYE|Geothlypis trichas|
|DERO|Desmodus rotundus|
|EPBR|Eptesicus brasiliensis|
|EPFU|Eptesicus fuscus|
|EPFUR|Eptesicus furinalis|
|EUGL|Eumops glaucinus|
|EUMA|Euderma maculatum|
|EUUN|Eumops underwoodi|
|GLSO|Glossophaga soricina|
|HAAM|Hemignathus virens|
|JABW|Horornis diphone|
|JAWE|Zosterops japonicus|
|LABL|Lasiurus blossevillii|
|LAEG|Lasiurus ega|
|LEGO|Spinus psaltria|
|LUWA|Oreothlypis luciae|
|MACR|_Carpodacus erythrinus_ used as substitute|
|MAWA|Setophaga magnolia|
|MELT|Garrulax canorus|
|MOMO|Molossus molossus|
|MOSI|Molossus sinaloae|
|MYNI|Myotis nigricans|
|MYRI|Myotis riparius|
|NOAL|Noctilio albiventris|
|NOLE|Noctilio leporinus|
|PHDI|Phyllostomus discolor|
|PHHA|Phyllostomus hastatus|
|PHST|Phylloderma stenops|
|PTGY|Pteronotus gymnonotus|
|PTPA|Pteronotus parnellii|
|PTPE|Pteronotus personatus|
|RBLE|Leiothrix lutea|
|RHBI|Rhogeessa bickhami|
|RHNA|Rhynchonycteris naso|
|RHTU|Rhogeessa tumida|
|SOSP|Melospiza melodia|
|STLI|Sturnira lilium|
|STPA|Sturnira parvidens|
|TOSA|Tonatia saurophila|
|TRCI|Trachops cirrhosus|
|VERD|Auriparus flaviceps|
|WCSP|Zonotrichia leucophrys|
|WIWA|Cardellina pusilla|
|YBCH|Icteria virens|
|YEWA|Setophaga petechia|

# QIIME formatting of host sequences
The above fasta sequences were converted into QIIME-formatted text and fasta files for import using the following script (imported concatenated fasta file starts as `tmp2.fa`)
```
seqtk seq -l0 tmp.fa > tmp2.fa
cat tmp2.fa | paste - - | awk -F "\t" '{print $2}' > tmp.seqs.txt
cat tmp2.fa | paste - - | awk -F "\t" '{print $1}' > tmp.header.txt
sed 's/>//g' tmp.header.txt | cut -f 1 -d ';' | awk '{ printf("%s_%s\n",$0,i++);next;} { print $0;}' > tax.accessions.txt
cut -d ';' -f 1 --complement tmp.header.txt | sed 's/tax=//g' > tmp.taxstring.txt


R
library(data.table)
setwd('/mnt/lustre/macmaneslab/devon/guano/BOLDdb')
df <- fread(file="tmp.taxstring.txt", header = FALSE, sep = ",", fill=TRUE)
df <- data.frame(lapply(df, function(x) {gsub("*.\\:", "", x)}))
df$V1 <- paste("k__", df$V1, sep = "")
df$V2 <- paste("p__", df$V2, sep = "")
df$V3 <- paste("c__", df$V3, sep = "")
df$V4 <- paste("o__", df$V4, sep = "")
df$V5 <- paste("f__", df$V5, sep = "")
df$V6 <- paste("g__", df$V6, sep = "")
df$V7 <- paste("s__", df$V7, sep = "")
df2 <- data.frame(paste(df$V1, df$V2, df$V3, df$V4, df$V5, df$V6, df$V7, sep = ";"))
write.table(df2, file = "/mnt/lustre/macmaneslab/devon/guano/BOLDdb/taxstring.txt", quote=FALSE, row.names = FALSE, col.names = FALSE)
q()
n

paste tax.accessions.txt taxstring.txt > qiime.host.txt
sed 's/^/>/g' tax.accessions.txt > tmp.tax.accessions.txt
paste tmp.tax.accessions.txt tmp.seqs.txt | tr '\t' '\n' > qiime.host.fa

rm tmp.seqs.txt tmp.header.txt tax.accessions.txt tmp.taxstring.txt tmp.tax.accessions.txt taxstring.txt
```

Files were then imported into QIIME-format with the following pair of functions:
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path qiime.host.fa \
  --output-path host_seqs.qza

qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path qiime.host.txt \
  --output-path host_taxonomy.qza
```

# References

Clare, E. L., Lim, B. K., Fenton, M. B., & Hebert, P. D. N. (2011). Neotropical Bats: Estimating Species Diversity with DNA Barcodes. PLoS ONE, 6(7), e22648. doi:10.1371/journal.pone.0022648 

Hebert,P.D., Stoeckle,M.Y., Zemlak,T.S. and Francis,C.M. (2004) Identification of Birds through DNA Barcodes. PLoS Biol. 2 (10), E312.

Kerr,K.C., Stoeckle,M.Y., Dove,C.J., Weigt,L.A., Francis,C.M. and Hebert,P.D. (2007) Comprehensive DNA barcode coverage of North American birds. Mol. Ecol. Notes 7(4), 535-543.

Kerr KC, Birks SM, Kalyakin MV, Red'kin YA, Koblik EA, Hebert PD. (2009) Filling the gap - COI barcode resolution in eastern Palearctic birds. Front Zool. Dec 9;6:29. doi: 10.1186/1742-9994-6-29.

Kerr,K.C. (2011) Searching for evidence of selection in avian DNA barcodes. Mol Ecol Resour 11 (6), 1045-1055.

Korstian JM, AM Hale, VJ Bennett, and DA Williams. 2016. Using DNA barcoding to improve bat carcass identification at wind farms in the United States. Conservation Genetics Resources 8:27-34.

Saitoh,T., Sugita,N., Someya,S., Iwami,Y., Kobayashi,S., Kamigaichi,H., Higuchi,A., Asai,S., Yamamoto,Y. and Nishiumi,I. (2015) DNA barcoding reveals 24 distinct lineages as cryptic bird species candidates in and around the Japanese Archipelago. Mol Ecol Resour 15 (1), 177-186.

Schindel DE, Stoeckle MY, Milensky C, Trizna M, Schmidt B, Gebhard C, Graves G. (2011) DNA barcodes of bird species in the national museum of natural history, smithsonian institution, USA.", Zookeys, Dec 8;(152):87-92.

Streicker DG, Turmelle AS, Vonhof MJ, Kuzmin IV, McCracken GF, Rupprecht CE. Host phylogeny constrains cross-species emergence and establishment of rabies virus in bats. Science. 2010 Aug 6;329(5992):676-9. doi: 10.1126/science.1188836.

Vonhof,M.J., Russell,A.L. and Miller-Butterworth,C.M. (07-08-2015) Range-Wide Genetic Analysis of Little Brown Bat (Myotis lucifugus) Populations: Estimating the Risk of Spread of White-Nose Syndrome. PLoS ONE 10:(7)E0128713


Walker FM, Williamson CHD, Sanchez DE, Sobek CJ, Chambers CL (2016) Species From Feces: Order-Wide Identification of Chiroptera From Guano and Other Non-Invasive Genetic Samples. PLOS ONE 11(9): e0162342. https://doi.org/10.1371/journal.pone.0162342
