COMPMUT RNAseq Analysis 7: Identifying lexA binding sites in SBW25
================
jpjh
compiled Feb 2021

## Identifying *lexA* regulated genes in *Pseudomonas fluorescens* SBW25

Identfied the following paper describing the SOS response in *P.
aeruginosa*.

[Cirz, Ryan T., Bryan M. O’Neill, Jennifer A. Hammond, Steven R. Head,
and Floyd E. Romesberg. 2006. “Defining the Pseudomonas Aeruginosa SOS
Response and Its Role in the Global Response to the Antibiotic
Ciprofloxacin.” Journal of Bacteriology 188
(20):7101–10.](https://www.ncbi.nlm.nih.gov/pubmed/17015649)

This paper identifies the consensus LexA box for PAO1 as

> `CTGTATAAATATACAG`

However only 6 of these residues are 100% conserved, instead the
sequence is:

> `CTG [GT][AC]TNN [AT][ACT][AGT][ACT][ACT] CAG`

Using ambiguous nucleotides the sequence is:

> `CTGKMTNNWHDHHCAG`

Use this to search the SBW25 chromosome. Re-purposing previously-written
scripts `sequenceMatches.py` and `getNextFeature.py`.

``` python
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqUtils

chr_file = "../rnaseq/ref/AM181176.gb"
lexA = "CTGKMTNNWHDHHCAG"

handle = open(chr_file, "rU")
reference = next(SeqIO.parse(handle, "gb"))
seq = str(reference.seq.upper())

q  = Seq(lexA.upper())
rq = q.reverse_complement()
matches_f = SeqUtils.nt_search(seq, str(q))
matches_r = SeqUtils.nt_search(seq, str(rq))

positions_forward = matches_f[1:len(matches_f)]
positions_reverse = matches_r[1:len(matches_r)]

table = []

for p in positions_forward:
   pos = int(p)
   for feature in reference.features:
       if feature.type == "CDS" and feature.strand == 1:
           if feature.location.start > pos:
               name = feature.qualifiers["locus_tag"][0]
               gene_start = feature.location.start
               difference = gene_start - pos
               strand = feature.strand
               table.append([pos, difference, name, strand])
               break

for q in positions_reverse:
   pos = int(q)
   for feature in reference.features:
       if feature.type == "CDS" and feature.strand == -1:
           if feature.location.end > pos:
               table.append(line)
               break
           name = feature.qualifiers["locus_tag"][0]
           gene_start = feature.location.end
           difference = pos - gene_start
           strand = feature.strand
           line = [pos, difference, name, strand]

output = open("../rnaseq/ref/LexA_sites.csv", "w+")

for line in table:
  output.write(",".join(str(x) for x in line))
  output.write("\n")
  
output.close()
```

Load the table in R.

``` r
lexA <- read.csv("../rnaseq/ref/LexA_sites.csv", sep=",", header=FALSE)
colnames(lexA) <- c("position","distance","locus_tag","strand")
lexA <- with(lexA, lexA[order(position),])
```

Identifies 36 putative lexA binding sites. However many of these are
some distance from the start of the following gene.

In Cirz et al., they identify 15 genes in the lexA regulon, all of which
are within 150 bp downstream of the lexA binding site. Applying this
criterion to our list gives a shortlist:

``` r
kable(lexA[lexA$distance<151,])
```

|     | position | distance | locus\_tag | strand |
|:----|---------:|---------:|:-----------|-------:|
| 1   |    54429 |       31 | PFLU\_0054 |      1 |
| 24  |   957428 |       39 | PFLU\_0848 |     -1 |
| 4   |  1008122 |       60 | PFLU\_0902 |      1 |
| 25  |  1303652 |      150 | PFLU\_1166 |     -1 |
| 6   |  1318057 |       57 | PFLU\_1189 |      1 |
| 7   |  1433371 |       70 | PFLU\_1295 |      1 |
| 26  |  1460562 |       61 | PFLU\_1323 |     -1 |
| 27  |  1709666 |       29 | PFLU\_1560 |     -1 |
| 10  |  1918050 |       42 | PFLU\_1750 |      1 |
| 28  |  3431180 |       19 | PFLU\_3139 |     -1 |
| 36  |  5421481 |       70 | PFLU\_4940 |     -1 |
| 21  |  5787859 |       40 | PFLU\_5271 |      1 |

------------------------------------------------------------------------

**[Back to index.](COMPMUT_index.md)**
