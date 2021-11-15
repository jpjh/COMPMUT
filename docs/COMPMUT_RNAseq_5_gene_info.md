COMPMUT RNAseq Analysis 5: Making Gene Info tables
================
jpjh
compiled Feb 2021

[Now published in PLoS
Biology](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001225):

Hall, J. P. J., Wright, R. C. T., Harrison, E., Muddiman, K. J., Jamie
Wood, A., Paterson, S., & Brockhurst, M. A. (2021). Plasmid fitness
costs are caused by specific genetic conflicts enabling resolution by
compensatory mutation. *PLoS Biology*, *19*(10), e3001225.
<https://doi.org/10.1371/journal.pbio.3001225>

**[Back to index.](COMPMUT_index.md)**

------------------------------------------------------------------------

## Getting gene\_info for chromosomal genes

Require a table with locus\_tag and details for each gene.

### Position information from the gff file

First: get details for each gene from the gff file. Use the gff file
extracted from the genbank using Python (described in
`COMPMUT_RNAseq_3_Rsubread.md`).

Load up functions for parsing GFF files, and use these to extract
relevant information from the chromosomal .gff file. From
[here](https://stat.ethz.ch/pipermail/bioconductor/2008-October/024669.html).

``` r
getAttributeField <- function (x, field, attrsep = ";") {
     s = strsplit(x, split = attrsep, fixed = TRUE)
     sapply(s, function(atts) {
         a = strsplit(atts, split = "=", fixed = TRUE)
         m = match(field, sapply(a, "[", 1))
         if (!is.na(m)) {
             rv = a[[m]][2]
         }
         else {
             rv = as.character(NA)
         }
         return(rv)
     })
}

gffRead <- function(gffFile, nrows = -1) {
     cat("Reading ", gffFile, ": ", sep="")
     gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
     header=FALSE, comment.char="#", nrows = nrows,
     colClasses=c("character", "character", "character", "integer",  
"integer",
     "character", "character", "character", "character"))
     colnames(gff) = c("seqname", "source", "feature", "start", "end",
             "score", "strand", "frame", "attributes")
     cat("found", nrow(gff), "rows with classes:",
         paste(sapply(gff, class), collapse=", "), "\n")
     stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
     return(gff)
}

gff <- gffRead("../rnaseq/ref/AM181176_python.gff")
```

    ## Reading ../rnaseq/ref/AM181176_python.gff: found 27897 rows with classes: character, character, character, integer, integer, character, character, character, character

``` r
gff$locus_tag <- getAttributeField(gff$attributes, "locus_tag")
gff$product <- getAttributeField(gff$attributes, "product")

gff_info <- gff[gff$feature=="CDS",c("start","end","strand","locus_tag","product")]
```

Remove NA locus\_tags from the gff\_info table.

``` r
gff_info <- gff_info[!(is.na(gff_info$locus_tag)),]
nrow(gff_info)
```

    ## [1] 6009

This gives 6009 CDS, consistent with Silby et al. (2009).

### GO terms from Uniprot

Collected details for each gene from Uniprot and saved as
`../rnaseq/ref/SBW25_uniprot.tab`. This tab file was run through a basic
AWK script to pull out only locus\_tags corresponding to SBW25
(i.e. PFLU\_), and piped out to `./ref/SBW25_uniprot.tsv`.

``` bash
cat ../rnaseq/ref/SBW25_uniprot.tab | awk -v FS="\t" '$7 ~ /^PFLU_/ {print $0}' \
  > ../rnaseq/ref/SBW25_uniprot.tsv
```

Load this file in R and examine.

``` r
uniprot_info <- read.table("../rnaseq/ref/SBW25_uniprot.tsv", sep="\t", header=FALSE, 
                           fill=TRUE, quote="", stringsAsFactors=FALSE)
colnames(uniprot_info) <- c("acc","entry_name","protein_name","gene_name","length","go","locus_tag")
nrow(uniprot_info)
```

    ## [1] 5914

There are fewer entries in the Uniprot file. See what entries from the
GFF are missing.

``` r
nrow(gff_info[!(gff_info$locus_tag %in% uniprot_info$locus_tag),])
```

    ## [1] 98

Almost 100 of these have no details from Uniprot (e.g. GO terms).

``` r
nrow(uniprot_info[!(uniprot_info$locus_tag %in% gff_info$locus_tag),])
```

    ## [1] 3

Three entries in the Uniprot database aren’t in the gff. Inspection
shows that this is due to a redundant numbering system that has one
entry for several identical sequences (e.g. insertion sequences).
Collect these sequences, remove them from the uniprot database, correct
them and re-add them as required.

``` r
duplicates <- uniprot_info[grep(" ",uniprot_info$locus_tag),]

uniprot_info <- uniprot_info[!(uniprot_info$acc %in% duplicates$acc),]

splitDuplicates <- function(x){
  genes <- strsplit(x$gene_name, split=" ")[[1]]
  outtab <- data.frame()
  for (i in 1:length(genes)){
    outtab <- rbind(outtab, x)
  }
  outtab$gene_name <- genes
  outtab$locus_tag <- genes
  return(outtab)
}

dups_expanded <- adply(duplicates, 1, function(x) splitDuplicates(x))

uniprot_info <- rbind(uniprot_info, dups_expanded)

nrow(uniprot_info)
```

    ## [1] 5920

Now just 89 are missing. Manual inspection of some of these on Uniprot
suggest that there aren’t Uniprot entries for these, with many of them
(87 pseudogenes). Pull all together and output.

``` r
gene_info <- merge(gff_info, uniprot_info, by="locus_tag", all=TRUE)
gene_info$length_na <- with(gene_info, (end-start)+1)

write.table(gene_info, file="../rnaseq/ref/AM181176_gene_info.tsv", quote=TRUE, row.names=FALSE, sep="\t")
```

## Getting gene info for plasmid genes

### pQBR57

For the plasmids, the .gff file were generated using Artemis.

``` r
gff_pq57 <- gffRead("../rnaseq/ref/LN713926.gff")
```

    ## Reading ../rnaseq/ref/LN713926.gff: found 426 rows with classes: character, character, character, integer, integer, character, character, character, character

``` r
gff_pq57$locus_tag <- getAttributeField(gff_pq57$attributes, "locus_tag")

gff_pq57_info <- gff_pq57[gff_pq57$feature=="CDS",c("start","end","strand","locus_tag")]
```

Remove NA tags as above.

``` r
gff_pq57_info <- gff_pq57_info[!(is.na(gff_pq57_info$locus_tag)),]
nrow(gff_pq57_info)
```

    ## [1] 426

Collect the uniprot info on pQBR57

``` bash
cat ../rnaseq/ref/SBW25_uniprot.tab | awk -v FS="\t" '$4 ~ /PQBR57_/ {print $0}' \
  > ../rnaseq/ref/pQBR57_uniprot.tsv
```

Load this file in R and examine.

``` r
uniprot_pq57_info <- read.table("../rnaseq/ref/pQBR57_uniprot.tsv", sep="\t", header=FALSE, 
                           fill=TRUE, quote="", stringsAsFactors=FALSE)
colnames(uniprot_pq57_info) <- c("acc","entry_name","protein_name","locus_tag","length","go")
nrow(uniprot_pq57_info)
```

    ## [1] 426

There are a handful of genes with no matching Uniprot entry…

``` r
nrow(gff_pq57_info[!(gff_pq57_info$locus_tag %in% uniprot_pq57_info$locus_tag),])
```

    ## [1] 8

These are mer operon genes and a recA homologue, and are probably due to
redundant entries in the same line of the table as described above. Sort
these out as before.

``` r
uniprot_pq57_info$gene_name <- uniprot_pq57_info$locus_tag
dups_pq57 <- uniprot_pq57_info[grep(" ",uniprot_pq57_info$locus_tag),]
uniprot_pq57_info <- uniprot_pq57_info[!(uniprot_pq57_info$acc %in% dups_pq57$acc),]

dups_pq57_expanded <- data.frame()
for(i in 1:nrow(dups_pq57)){
  dups_pq57_expanded <- rbind(dups_pq57_expanded, splitDuplicates(dups_pq57[i,]))
}

dups_pq57_expanded <- dups_pq57_expanded[grep("^PQBR57_", dups_pq57_expanded$locus_tag),]

uniprot_pq57_info <- rbind(uniprot_pq57_info,dups_pq57_expanded)
nrow(uniprot_pq57_info)
```

    ## [1] 426

Pull together and output.

``` r
pq57_gene_info <- merge(gff_pq57_info, uniprot_pq57_info, by="locus_tag", all=TRUE)

write.table(pq57_gene_info, file="../rnaseq/ref/LN713926_gene_info.tsv", quote=TRUE, row.names=FALSE, sep="\t")
```

### pQBR103

``` r
gff_pq103 <- gffRead("../rnaseq/ref/AM235768.gff")
```

    ## Reading ../rnaseq/ref/AM235768.gff: found 1052 rows with classes: character, character, character, integer, integer, character, character, character, character

``` r
gff_pq103$locus_tag <- getAttributeField(gff_pq103$attributes, "locus_tag")
gff_pq103_info <- gff_pq103[gff_pq103$feature=="CDS",c("start","end","strand","locus_tag")]
gff_pq103_info <- gff_pq103_info[!(is.na(gff_pq103_info$locus_tag)),]
nrow(gff_pq103_info)
```

    ## [1] 479

Collect the uniprot info on pQBR103.

``` bash
cat ../rnaseq/ref/SBW25_uniprot.tab | awk -v FS="\t" '$7 ~ /pQBR0/ {print $0}' \
  > ../rnaseq/ref/pQBR103_uniprot.tsv
```

Load this file in R and examine.

``` r
uniprot_pq103_info <- read.table("../rnaseq/ref/pQBR103_uniprot.tsv", sep="\t", header=FALSE, 
                           fill=TRUE, quote="", stringsAsFactors=FALSE)
colnames(uniprot_pq103_info) <- c("acc","entry_name","protein_name","locus_tag","length","go","gene_name")
nrow(uniprot_pq103_info)
```

    ## [1] 474

There are a handful of genes with no matching Uniprot entry…

``` r
nrow(gff_pq103_info[!(gff_pq103_info$locus_tag %in% uniprot_pq103_info$locus_tag),])
```

    ## [1] 31

Various genes, arising due to double entries in the Uniprot database
again.

``` r
uniprot_pq103_info$gene_name <- uniprot_pq103_info$locus_tag
dups_pq103 <- uniprot_pq103_info[grep(" ",uniprot_pq103_info$locus_tag),]
uniprot_pq103_info <- uniprot_pq103_info[!(uniprot_pq103_info$acc %in% dups_pq103$acc),]

dups_pq103_expanded <- data.frame()

for(i in 1:nrow(dups_pq103)){
  dups_pq103_expanded <- rbind(dups_pq103_expanded, splitDuplicates(dups_pq103[i,]))
}

dups_pq103_expanded <- dups_pq103_expanded[grep("^pQBR0", dups_pq103_expanded$locus_tag),]

uniprot_pq103_info <- rbind(uniprot_pq103_info,dups_pq103_expanded)
nrow(uniprot_pq103_info)
```

    ## [1] 474

Still missing 5 genes.

``` r
gff_pq103_info[!(gff_pq103_info$locus_tag %in% uniprot_pq103_info$locus_tag),]
```

    ##      start    end strand locus_tag
    ## 37   14242  14496      +  pQBR0016
    ## 262 129713 130228      -  pQBR0120
    ## 673 287675 288058      +  pQBR0308
    ## 936 377569 377868      +  pQBR0429
    ## 937 377869 377940      +  pQBR0429

One of these, pQBR0429 is a pseudogene, and it has two entries for each
frame that it is in.

Pull together and output. Correct the ‘length’ column so all entries
have a gene length.

``` r
pq103_gene_info <- merge(gff_pq103_info, uniprot_pq103_info, by="locus_tag", all=TRUE)
pq103_gene_info$length <- with(pq103_gene_info, (end-start+1))

# remove duplicate entries, i.e. pQBR429

pq103_gene_info <- pq103_gene_info[!duplicated(pq103_gene_info$locus_tag),]

write.table(pq103_gene_info, file="../rnaseq/ref/AM235768_gene_info.tsv", quote=TRUE, row.names=FALSE, sep="\t")
```

------------------------------------------------------------------------

**[Back to index.](COMPMUT_index.md)**
