COMPMUT RNAseq Analysis 2: Alignment of reads to the reference genome
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

## Aligning reads to the reference

Using `HISAT2` for alignment, as at the time of analysis `bwa` and
`bowtie2` do not obviously allow strand-specific alignment.

Build the references. Each treatment will be mapped onto the
corresponding reference sequence(s).

The EBI accession numbers for the different sequences are as follows:

-   SBW25 = AM181176
-   PQBR57 = LN713926
-   PQBR103 = AM235768

``` bash
mkdir ../rnaseq/ref

curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=AM181176&amp;rettype=fasta" \
   > ../rnaseq/ref/AM181176.fasta
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=LN713926&amp;rettype=fasta" \
   > ../rnaseq/ref/LN713926.fasta
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&amp;id=AM235768&amp;rettype=fasta" \
   > ../rnaseq/ref/AM235768.fasta

~/programs/hisat2-2.1.0/hisat2-build \
   -f ../rnaseq/ref/AM181176.fasta \
   ../rnaseq/ref/SBW25

~/programs/hisat2-2.1.0/hisat2-build \
   -f ../rnaseq/ref/AM181176.fasta,../rnaseq/ref/LN713926.fasta \
   ../rnaseq/ref/SBW25_pQBR57

~/programs/hisat2-2.1.0/hisat2-build \
   -f ../rnaseq/ref/AM181176.fasta,../rnaseq/ref/AM235768.fasta \
   ../rnaseq/ref/SBW25_pQBR103
```

Note that for the [NEB Ultra Directional
kit](https://international.neb.com/products/e7420-nebnext-ultra-directional-rna-library-prep-kit-for-illumina)
that we used for this protocol, we use `rna-strandness RF`.

For bacteria, use no-spliced-alignment. Exploratory studies showed that
including splicing resulted in abundant spurious splice events, and
[since splicing in bacteria is thought to be
rare](https://www.nature.com/scitable/topicpage/rna-splicing-introns-exons-and-spliceosome-12375/)
(but see [here](http://genesdev.cshlp.org/content/12/9/1243.full.html))
it’s unlikely to make a big contribution to overall transcription
patterns.

``` r
kable(read.table("../rnaseq/TreatmentTable.csv", header=TRUE, sep=","))
```

| trt | plasmid | amelioration        |
|:----|:--------|:--------------------|
| R01 | none    | wt                  |
| R02 | none    | dPFLU4242           |
| R03 | pQBR57  | wt                  |
| R04 | pQBR57  | dPFLU4242           |
| R05 | pQBR103 | wt                  |
| R06 | pQBR103 | dPFLU4242           |
| R07 | pQBR57  | PQBR57\_0059\_V100A |
| R08 | none    | dgacS               |
| R09 | pQBR57  | dgacS               |
| R10 | pQBR103 | dgacS               |

Sample names ending 1, 2, 8 are mapped against SBW25.

Sample names ending 3, 4, 7, 9 are mapped against SBW25\_pQBR57.

Sample names ending 5, 6, 0 are mapped against SBW25\_pQBR103.

``` bash
mkdir ../rnaseq/hisat_summary

find ../rnaseq/reads/Trimmed | grep "L001_R1_001.fastq.gz$" | sort > ../rnaseq/L1R1reads.txt
find ../rnaseq/reads/Trimmed | grep "L001_R2_001.fastq.gz$" | sort > ../rnaseq/L1R2reads.txt
find ../rnaseq/reads/Trimmed | grep "L002_R1_001.fastq.gz$" | sort > ../rnaseq/L2R1reads.txt
find ../rnaseq/reads/Trimmed | grep "L002_R2_001.fastq.gz$" | sort > ../rnaseq/L2R2reads.txt
ls -1 ../rnaseq/reads/Trimmed | sort | tr -d "/" > ../rnaseq/samples.txt

paste ../rnaseq/samples.txt \
  ../rnaseq/L1R1reads.txt \
  ../rnaseq/L1R2reads.txt \
  ../rnaseq/L2R1reads.txt \
  ../rnaseq/L2R2reads.txt > ../rnaseq/sampreads.txt

grep -E "Sample_[0-9]+-[ABC]R[01][128]" ../rnaseq/sampreads.txt \
  | while read sample L1R1 L1R2 L2R1 L2R2
do
~/programs/hisat2-2.1.0/hisat2 \
       --phred33 \
       -q -p 8 \
       --no-spliced-alignment \
       --summary-file ../rnaseq/hisat_summary/${sample}_hisat.txt \
       --new-summary \
       --rna-strandness RF \
       -x ../rnaseq/ref/SBW25 \
       -1 $L1R1,$L2R1 -2 $L1R2,$L2R2 \
    | samtools view -q10 -bSh \
    > ../rnaseq/aligned.bam
samtools view -b ../rnaseq/aligned.bam \
    | samtools sort -@ 8 -m 1G > ../rnaseq/bams/${sample}.bam
samtools index ../rnaseq/bams/${sample}.bam
done

grep -E "Sample_[0-9]+-[ABC]R[01][3479]" ../rnaseq/sampreads.txt \
  | while read sample L1R1 L1R2 L2R1 L2R2
do
~/programs/hisat2-2.1.0/hisat2 \
       --phred33 \
       -q -p 8 \
       --no-spliced-alignment \
       --summary-file ../rnaseq/hisat_summary/${sample}_hisat.txt \
       --new-summary \
       --rna-strandness RF \
       -x ../rnaseq/ref/SBW25_pQBR57 \
       -1 $L1R1,$L2R1 -2 $L1R2,$L2R2 \
    | samtools view -q10 -bSh \
    > ../rnaseq/aligned.bam
samtools view -b ../rnaseq/aligned.bam \
    | samtools sort -@ 8 -m 1G > ../rnaseq/bams/${sample}.bam
samtools index ../rnaseq/bams/${sample}.bam
done

grep -E "Sample_[0-9]+-[ABC]R[01][560]" ../rnaseq/sampreads.txt  \
  | while read sample L1R1 L1R2 L2R1 L2R2
do
~/programs/hisat2-2.1.0/hisat2 \
       --phred33 \
       -q -p 8 \
       --no-spliced-alignment \
       --summary-file ../rnaseq/hisat_summary/${sample}_hisat.txt \
       --new-summary \
       --rna-strandness RF \
       -x ../rnaseq/ref/SBW25_pQBR103 \
       -1 $L1R1,$L2R1 -2 $L1R2,$L2R2 \
    | samtools view -q10 -bSh \
    > ../rnaseq/aligned.bam
samtools view -b ../rnaseq/aligned.bam \
    | samtools sort -@ 8 -m 1G > ../rnaseq/bams/${sample}.bam
samtools index ../rnaseq/bams/${sample}.bam
done
```

### Viewing coverage by strand

Note that by selecting the ‘Colour By / RNAseq Strand Specific Tag (XS)’
and using the ‘Strand Coverage’ option, it’s possible to view several
coverage plots in Artemis very conveniently. They can be coloured using
the ‘Views / Coverage Options / Configure Lines…’ function. However,
note that these aren’t scaled by total read counts, so it’s difficult to
compare between samples.

Extracting the data from the .bam file can be done using the `XS:A` tag
added by [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml).
Alternatively, since the first read in each pair is opposite sense to
the RNA, and the second read is sense to the RNA, this information can
be extracted using the flag specification as detailed [in this Biostars
post](https://www.biostars.org/p/92935/): *forward reads* are selected
using `-f 128 -F 16` (second in pair, not read reverse strand) and
`-f 80` for reverse (read reverse strand, first in pair); *reverse
reads* are selected using `-f 144` (second in pair, reverse) and
`-f 64 -F16` (first in pair, not read reverse).

Modified the Biostars script to give `strandsplit_bam_cov.sh`, which
also extracts coverage from each strand using genomeCoverageBed.

    # strandsplit_bam_cov.sh
    # Script for getting sense and antisense reads from a bam file
    # Assuming that library was made with reverse-directional prep (e.g. UltraDirectional)
    # Based on script found here https://www.biostars.org/p/92935/

    set -uex

    # Get the bam file from the command line
    DATA=$1
    OUTPUT=$2

    # Forward strand.
    # 1. alignments of the second in pair if they map to the forward strand
    # 2. alignments of the first in pair if they map to the reverse strand

    samtools view -b -f 128 -F 16 $DATA > fwd1.temp.bam
    samtools index fwd1.temp.bam

    samtools view -b -f 80 $DATA > fwd2.temp.bam
    samtools index fwd2.temp.bam

    # Combine alignments that originate on the forward strand.
    samtools merge -f fwd.temp.bam fwd1.temp.bam fwd2.temp.bam
    samtools index fwd.temp.bam

    # Reverse strand
    # 1. alignments of the second in pair if they map to the reverse strand
    # 2. alignments of the first in pair if they map to the forward strand
    samtools view -b -f 144 $DATA > rev1.temp.bam
    samtools index rev1.temp.bam

    samtools view -b -f 64 -F 16 $DATA > rev2.temp.bam
    samtools index rev2.temp.bam

    # Combine alignments that originate on the reverse strand.
    samtools merge -f rev.temp.bam rev1.temp.bam rev2.temp.bam
    samtools index rev.temp.bam

    # Run genomeCoverageBed on the forward and reverse bams

    genomeCoverageBed \
        -ibam fwd.temp.bam \
        -bga > ${OUTPUT}_f.cov.txt

    genomeCoverageBed \
        -ibam rev.temp.bam \
        -bga > ${OUTPUT}_r.cov.txt

Run this script on the `.bam` files.

``` bash
cat ../rnaseq/sampreads.txt | while read sample L1R1 L1R2 L2R1 L2R2
do
  bash ../functions/strandsplit_bam_cov.sh \
    ../rnaseq/bams/${sample}.bam \
    ../rnaseq/cov/${sample}
done
```

Note that Rsubread can extract the strandedness information
automatically, so this approach is only needed for plotting per-strand
coverage, i.e. the output from `genomeCoverageBed`.

[Now proceed to using Rsubread to calculate
counts-per-feature.](COMPMUT_RNAseq_3_Rsubread.md)

------------------------------------------------------------------------

**[Back to index.](COMPMUT_index.md)**
