COMPMUT RNAseq Analysis 3: Using Rsubread’s featureCounts to count reads
mapped to features
================
jpjh
compiled Feb 2021

## Analysis of alignments using Rsubread

First, need to generate a .gtf file for analysis with Rsubread.

``` bash
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=AM181176&amp;rettype=gb" \
   > ../rnaseq/ref/AM181176.gb
   
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=LN713926&amp;rettype=gb" \
   > ../rnaseq/ref/LN713926.gb
   
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=AM235768&amp;rettype=gb" \
   > ../rnaseq/ref/AM235768.gb
```

Convert these files to .gff using Artemis, then run the following bash
code to get a .gtf file.

``` bash
cat ../rnaseq/ref/AM181176.gff | awk '$3=="CDS" {print $0}' \
    | sed 's/AM181176/AM181176.4/g' \
    | awk -v FS=";" '{print $1}' | sed 's/ID=//g' \
    | awk -v OFS="\t" \
          '{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \"" $9 "\"; transcript_id \"" $9 "\""}' \
          > ../rnaseq/ref/AM181176.gtf
                    
cat ../rnaseq/ref/LN713926.gff | awk '$3=="CDS" {print $0}' \
    | sed 's/LN713926/LN713926.1/g' \
    | awk -v FS=";" '{print $1}' | sed 's/ID=//g' \
    | awk -v OFS="\t" \
          '{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \"" $9 "\"; transcript_id \"" $9 "\""}' \
          > ../rnaseq/ref/LN713926.gtf
          
cat ../rnaseq/ref/AM235768.gff | awk '$3=="CDS" {print $0}' \
    | sed 's/AM235768/AM235768.1/g' \
    | awk -v FS=";" '{print $1}' | sed 's/ID=//g' \
    | awk -v OFS="\t" \
          '{print $1,$2,$3,$4,$5,$6,$7,$8,"gene_id \"" $9 "\"; transcript_id \"" $9 "\""}' \
          > ../rnaseq/ref/AM235768.gtf
          
cat ../rnaseq/ref/AM181176.gtf \
  ../rnaseq/ref/LN713926.gtf \
  ../rnaseq/ref/AM235768.gtf > ../rnaseq/ref/SBW25_pQBR57_pQBR103.gtf

mkdir ../rnaseq/rsubread
```

Install required libraries for Rsubread.

``` r
BiocManager::install("Rsubread")
```

Load up the required libraries.

``` r
library("Rsubread")
sample_ids <- dir(path="../rnaseq/bams", full.names=TRUE, pattern="*.bam$")
```

Run `featureCounts`.

Using `strandSpecific = 2` since the library prep is reverse stranded
(exploratory analysis using Sample\_1 showed that this was the correct
option). Specifying paired-end reads but not requiring both ends to map
(I figure that since some genes are transcribed as operons it is
possible that the mate for one read might map to the adjacent feature).

``` r
counts <- featureCounts(sample_ids, annot.ext = "../rnaseq/ref/SBW25_pQBR57_pQBR103.gtf",
                        isGTFAnnotationFile = TRUE,
                        GTF.featureType = "CDS",
                        GTF.attrType = "gene_id",
                        allowMultiOverlap = TRUE,
                        isPairedEnd = TRUE,
                        strandSpecific = 2,
                        nthreads=1)

saveRDS(counts, file="../rnaseq/rsubread/counts_v1.rds")

write.table(counts$counts,
            file="../rnaseq/rsubread/CountsMatrix_v1.csv", 
            quote=FALSE, row.names=TRUE, sep=",")
write.table(counts$stat,
            file="../rnaseq/rsubread/CountsMatrix_v1_stats.csv", 
            quote=FALSE, row.names=TRUE, sep=",")
```

This produces counts matrices for each of the mappings. Check the
statistics for sanity.

``` r
stats <- data.frame(Status=counts$stat$Status, 
             total_mapped=rowSums(counts$stat[,-1]))
kable(stats)
```

| Status                          | total\_mapped |
|:--------------------------------|--------------:|
| Assigned                        |     445350255 |
| Unassigned\_Unmapped            |             0 |
| Unassigned\_MappingQuality      |             0 |
| Unassigned\_Chimera             |             0 |
| Unassigned\_FragmentLength      |             0 |
| Unassigned\_Duplicate           |             0 |
| Unassigned\_MultiMapping        |             0 |
| Unassigned\_Secondary           |             0 |
| Unassigned\_Nonjunction         |             0 |
| Unassigned\_NoFeatures          |     113393386 |
| Unassigned\_Overlapping\_Length |             0 |
| Unassigned\_Ambiguity           |             0 |

Reads are either assigned, or they are unassigned because there wasn’t a
feature in the .gtf.

``` r
with(stats, stats[Status=="Assigned","total_mapped"] / sum(stats[,"total_mapped"]))
```

    ## [1] 0.7970565

About 80% of reads are associated with features, with &gt;20% not
associated with any features. A considerable proportion of these will be
the reads mapping to ssrA (identified during the QC step), since
`misc_RNA` features weren’t transferred to the .gff file. Note that
Artemis doesn’t output `misc_RNA`, `regulatory`, `sig_peptide`,
`source`, `protein_bind`, `misc_feature`, `stem_loop`. Of these, only
`misc_RNA` identifies a gene (the others identify UTR of transcripts,
regulatory regions, or domains).

Reattempt assignment with a .gtf file with `misc_RNA` included.

``` bash
grep "misc_RNA" -A 1 ../rnaseq/ref/AM181176.gb
```

    ##      misc_RNA        557909..558014
    ##                      /note="TPP riboswitch (THI element)"
    ## --
    ##      misc_RNA        complement(1019665..1019806)
    ##                      /note="yybP-ykoY element"
    ## --
    ##      misc_RNA        1040892..1041184
    ##                      /note="Bacterial Ribonuclease P class A"
    ## --
    ##      misc_RNA        1299177..1299340
    ##                      /note="PrrB/RsmZ RNA family"
    ## --
    ##      misc_RNA        complement(2943100..2943306)
    ##                      /note="Cobalamin riboswitch"
    ## --
    ##      misc_RNA        2944182..2944381
    ##                      /note="Cobalamin riboswitch"
    ## --
    ##      misc_RNA        4296931..4297080
    ##                      /note="PrrF RNA"
    ## --
    ##      misc_RNA        4663279..4663474
    ##                      /note="Cobalamin riboswitch"
    ## --
    ##      misc_RNA        complement(4959120..4959373)
    ##                      /note="Cobalamin riboswitch"
    ## --
    ##      misc_RNA        complement(5051832..5051931)
    ##                      /note="Bacterial signal recognition particle RNA"
    ## --
    ##      misc_RNA        5226269..5226438
    ##                      /note="FMN riboswitch (RFN element)"
    ## --
    ##      misc_RNA        5735279..5735433
    ##                      /note="PrrF RNA"
    ## --
    ##      misc_RNA        5798308..5798701
    ##                      /note="Bacterial tmRNA, 10Sa RNA or SsrA"
    ## --
    ##      misc_RNA        complement(6141773..6141890)
    ##                      /note="RsmY RNA family"
    ## --
    ##      misc_RNA        6438562..6438739
    ##                      /note="6S / SsrS RNA"

Of these 15 hits, the riboswitches are cis-acting elements in mRNA so
are unlikely to vary greatly from the corresponding transcripts. This
leaves:

-   1040892..1041184 “Bacterial Ribonuclease P class A”
-   1299177..1299340 “PrrB/RsmZ RNA family”
-   4296931..4297080 “PrrF RNA”
-   complement(5051832..5051931) “Bacterial signal recognition particle
    RNA”
-   5735279..5735433 “PrrF RNA”
-   5798308..5798701 “Bacterial tmRNA, 10Sa RNA or SsrA”
-   complement(6141773..6141890) “RsmY RNA family”
-   6438562..6438739 “6S / SsrS RNA”

Make these into a table and append to
`../rnaseq/ref/SBW25_pQBR57_pQBR103.gtf`. Used the following general
code to create `AM181176_python.gff` from `AM181176.gb`

``` python
from BCBio import GFF
from Bio import SeqIO

in_file = "../rnaseq/ref/AM181176.gb"
out_file = "../rnaseq/ref/AM181176_python.gff"
in_handle = open(in_file)
out_handle = open(out_file, "w")

GFF.write(SeqIO.parse(in_handle, "genbank"), out_handle)

in_handle.close()
out_handle.close()
```

Run bash script to pull out the important lines…

``` bash
awk '$3=="misc_RNA" {print $0}' ../rnaseq/ref/AM181176_python.gff \
  | grep -v "riboswitch" \
  | awk -v FS=";" '{print $1}' \
    > ../rnaseq/ref/AM181176_misc_RNA.gff
```

…but I ended up having to edit it manually anyway to get the final
column correct. The edited file in
`../rnaseq/ref/AM181176_misc_RNA.gtf`. Compile this with the rest of the
CDS, and replace ‘CDS’ and ‘misc\_RNA’ with ‘exon’ to ensure these
features can be analysed together.

``` bash
cat ../rnaseq/ref/SBW25_pQBR57_pQBR103.gtf ../rnaseq/ref/AM181176_misc_RNA.gtf \
  | awk -v FS=" " \
        -v OFS="\t" \
        '{print $1, $2, "exon", $4, $5, $6, $7, $8, $9 " " $10 " " $11 " " $12}' \
  > ../rnaseq/ref/SBW25_pQBR57_pQBR103_RNA.gtf
```

Run through all samples.

``` r
counts_exon <- featureCounts(sample_ids, annot.ext = "../rnaseq/ref/SBW25_pQBR57_pQBR103_RNA.gtf",
                        isGTFAnnotationFile = TRUE,
                        GTF.featureType = "exon",
                        GTF.attrType = "gene_id",
                        allowMultiOverlap = TRUE,
                        isPairedEnd = TRUE,
                        strandSpecific = 2,
                        nthreads=1)

saveRDS(counts_exon, file="../rnaseq/Rsubread/counts.rds")

write.table(counts_exon$counts,
            file="../rnaseq/rsubread/CountsMatrix.csv", 
            quote=FALSE, row.names=TRUE, sep=",")
write.table(counts_exon$stat,
            file="../rnaseq/rsubread/CountsMatrix_stats.csv", 
            quote=FALSE, row.names=TRUE, sep=",")
```

See how these results compare with the previous ones.

``` r
stats_v2 <- data.frame(Status=counts_exon$stat$Status, 
                       total_mapped=rowSums(counts_exon$stat[,-1]))
kable(stats_v2)
```

| Status                          | total\_mapped |
|:--------------------------------|--------------:|
| Assigned                        |     531711574 |
| Unassigned\_Unmapped            |             0 |
| Unassigned\_MappingQuality      |             0 |
| Unassigned\_Chimera             |             0 |
| Unassigned\_FragmentLength      |             0 |
| Unassigned\_Duplicate           |             0 |
| Unassigned\_MultiMapping        |             0 |
| Unassigned\_Secondary           |             0 |
| Unassigned\_Nonjunction         |             0 |
| Unassigned\_NoFeatures          |      27032067 |
| Unassigned\_Overlapping\_Length |             0 |
| Unassigned\_Ambiguity           |             0 |

Calcalate proportion.

``` r
with(stats_v2, stats_v2[Status=="Assigned","total_mapped"] / sum(stats_v2[,"total_mapped"]))
```

    ## [1] 0.9516199

Now &gt;95% of reads are associated with a feature.

Check each bamfile for proportion assigned.

``` r
stats_bamfile <- counts_exon$stat %>% 
  filter(Status %in% c("Assigned", "Unassigned_NoFeatures"))
colnames(stats_bamfile) <- gsub("X.Users.jpjh.Dropbox.WORK.PROJECTS.18.09_COMPMUT.COMPMUT_RNAseq.bams.",
                                "", colnames(stats_bamfile))
stats_bamfile_t <- data.frame(t(stats_bamfile[,-1]))
colnames(stats_bamfile_t) <- c("Assigned","Unassigned_NoFeatures")
stats_bamfile_t$proportion <- with(stats_bamfile_t, Assigned/(Assigned + Unassigned_NoFeatures))
summary(stats_bamfile_t$proportion)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##  0.9386  0.9476  0.9514  0.9510  0.9533  0.9650

All samples have between 93.5% and 96.5% of reads assigned now that the
non-coding RNAs are included.

Can now [proceed to edgeR to perform differential
expression.](COMPMUT_RNAseq_4_edgeR.md)

------------------------------------------------------------------------

**[Back to index.](COMPMUT_index.md)**
