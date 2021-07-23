COMPMUT RNAseq Analysis 6c: pQBR103 genes — plots, figures, and analysis
================
jpjh
compiled Feb 2021, edited Jul 2021

## Analysis of edgeR tables outline

Recall:

This data was collected in
[`COMPMUT_RNAseq_4_edgeR.Rmd`](COMPMUT_RNAseq_4_edgeR.md).

This analysis will pull out differentially expressed genes under the
different treatments.

1.  Chromosomal genes differentially expressed relative to the ancestor.

-   Effect of plasmid acquisition
-   Effect of ameliorations without plasmid
-   Effect of ameliorations on effect of plasmid

1.  Chromosomal genes differentially expressed relative to the
    unameliorated plasmid.

-   Are the genes which are upregulated by the plasmid and then not
    significantly different from the ancestor significantly
    downregulated relative to the plasmid-carrying strain?

1.  Plasmid genes differentially expressed relative to the unameliorated
    plasmid.

-   Comparisons of gene expression from plasmid and from chromosome.

See [`COMPMUT_RNAseq_6a_analysis`](COMPMUT_RNAseq_6a_analysis.md) for
points 1 and 2. This document addresses point 3 for pQBR103. See
[`COMPMUT_RNAseq_6b_analysis`](COMPMUT_RNAseq_6b_analysis.md) for pQBR57
genes.

### 3b. Plasmid genes (pQBR103) differentially expressed relative to the unameliorated plasmid.

Read in tables of plasmid differential gene expression.

``` r
full_table_pq103 <- read.csv("../data/COMPMUT_RNAseq_4_de_table_pQ103.csv", sep=",", header=TRUE)
full_table_pq103 <- mutate(full_table_pq103, gene = rownames(full_table_pq103))
dim(full_table_pq103)
```

    ## [1] 478  11

Generate volcano plots for each treatment. Reshape data into long
format.

``` r
source("../functions/makeLong.R")

pq103 <- makeLong(full_table_pq103)
head(pq103)
```

    ##       gene          contrast        FDR   logCPM       logFC plasmid
    ## 1 pQBR0001     pQBR103.dgacS 0.04671859 14.18554  0.34785103 pQBR103
    ## 2 pQBR0001 pQBR103.dPFLU4242 0.75485616 14.18554  0.06234590 pQBR103
    ## 3 pQBR0002     pQBR103.dgacS 0.04571335 14.28672  0.35954675 pQBR103
    ## 4 pQBR0002 pQBR103.dPFLU4242 0.50259910 14.28672  0.12461086 pQBR103
    ## 5 pQBR0003     pQBR103.dgacS 0.10707498 10.83910 -0.40848857 pQBR103
    ## 6 pQBR0003 pQBR103.dPFLU4242 0.75783900 10.83910 -0.09034341 pQBR103
    ##   amelioration
    ## 1        dgacS
    ## 2    dPFLU4242
    ## 3        dgacS
    ## 4    dPFLU4242
    ## 5        dgacS
    ## 6    dPFLU4242

Generate volcano plots.

``` r
ggplot(data=pq103, aes(x=logFC, y=(-log10(FDR)))) +
  ggtitle("pQBR103 sense transcripts") +
  geom_vline(xintercept=0, size=0.1) +
  geom_vline(xintercept=c(-1,1), size=0.1, linetype="dashed") +
  geom_hline(yintercept=0, size=0.1) + 
  geom_hline(yintercept=-log10(0.05), size=0.1, linetype="dashed") + 
  geom_point(alpha=0.2, shape=16) +
  facet_grid(.~amelioration) 
```

![](COMPMUT_RNAseq_6c_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

#### Investigate effects of different ameliorations on expression

``` r
pq103 %>% group_by(amelioration) %>%
  filter(FDR<0.05) %>%
  summarise(de_all = n(),
            up_all = sum(logFC>0),
            down_all = sum(logFC<0),
            up_2x = sum(logFC>1),
            down_2x = sum(logFC<(-1))) %>% kable()
```

| amelioration | de\_all | up\_all | down\_all | up\_2x | down\_2x |
|:-------------|--------:|--------:|----------:|-------:|---------:|
| dPFLU4242    |      85 |      67 |        18 |     21 |        0 |
| dgacS        |      89 |      42 |        47 |      4 |        3 |

##### Plot effects of different ameliorations against one another

Load up gene info, remove any duplicates, and create a ‘gene’ column for
merging.

``` r
pq103_gene_info <- read.table("../rnaseq/ref/AM235768_gene_info.tsv", header=TRUE, sep="\t",
                      stringsAsFactors = FALSE)

pq103_gene_info <- pq103_gene_info[!(duplicated(pq103_gene_info)),] %>% mutate(gene = gene_name)
```

Load up information on the regions or regulons to which these genes
belong. From [Hall et al. 2015](dx.doi.org/10.1111/1462-2920.12901).
Define how they should be coloured in the figures. Merge into full table
and into the table of gene info.

``` r
pregions <- read.csv("../ext/pQBR_regions.csv", header=TRUE, sep=",") %>%
  mutate(region = ifelse(region=="", "N/A", region),
         gene = locus_tag)

regions_cols <- data.frame(
  region = c("tra",       "par",        "che",        "sam",    
             "pil",   "Tn5042", "uvr",    "rep",  "N/A"),
  colour = c("steelblue1","darkorange1","chartreuse2","yellow3",
             "purple","red3",   "hotpink","grey", "grey20")
)

pregions <- pregions %>% left_join(regions_cols, by="region")

pq103_gene_info <- left_join(pq103_gene_info, pregions, by="gene")

full_table_pq103 <- left_join(full_table_pq103, pregions, by="gene")
```

###### PFLU4242 and dgacS

First: the contrasting plasmid and chromosomal mutations.

``` r
pq103_de_am <- filter(pq103, FDR<0.05) 
pq103_de_am_list <- pq103_de_am %>% 
  select(gene) %>% pull() %>% unique()

tab_pq103_am <- full_table_pq103 %>% filter(gene %in% pq103_de_am_list) %>%
  select(matches("pQBR103\\.(dgacS|dPFLU4242)\\."), region, colour, gene)

tab_pq103_am <- tab_pq103_am %>% 
  mutate(sig = ifelse(pQBR103.dgacS.FDR<0.05 & pQBR103.dPFLU4242.FDR<0.05,"both",
                      ifelse(pQBR103.dgacS.FDR<0.05 & pQBR103.dPFLU4242.FDR>0.05,"dgacS only",
                             ifelse(pQBR103.dgacS.FDR>0.05 & pQBR103.dPFLU4242.FDR<0.05,"dPFLU4242 only",
                                    "neither"))))

leg_title <- "FDR <0.05"

lines <- data.frame(intercept= c(-1,0,1),
                    linetype = c("2","1","2"))

(p_pq103_am <- ggplot(data=tab_pq103_am,
       aes(x=pQBR103.dgacS.logFC, y=pQBR103.dPFLU4242.logFC)) +
  geom_vline(data=lines, aes(xintercept=intercept, linetype=linetype), size=0.2, show.legend=FALSE) +
  geom_hline(data=lines, aes(yintercept=intercept, linetype=linetype), size=0.2, show.legend=FALSE) +
  geom_point(aes(shape=sig, size=sig, colour=region), alpha=0.8) +
  scale_size_manual(values=c(1.5,0.8,0.8,0.8)) +
  scale_shape_manual(values=c(16,0,4,16)) +
  scale_colour_manual(values=regions_cols$colour, breaks=regions_cols$region) +
  coord_fixed() +
  labs(x="pQBR103 with ∆gacS mutation log fold-change",
       y="pQBR103 with ∆PFLU4242 log fold-change",
       alpha=leg_title, size=leg_title, shape=leg_title) +
  theme(legend.position="right"))
```

![](COMPMUT_RNAseq_6c_analysis_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

Plot different significances on different plots.

``` r
p_pq103_am + facet_wrap(~sig)
```

![](COMPMUT_RNAseq_6c_analysis_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

The pattern for pQBR103 is very different from that of pQBR57, with a
large number of upregulated genes in both compensated variants. These
include the transfer genes, and the pilus genes which are downregulated
in pQBR57.

Genes downregulated by both are:

``` r
pq103_am_both <- tab_pq103_am %>% filter(sig=="both") %>% select(gene) %>% pull() %>% unique()

pq103_gene_info %>% filter(gene_name %in% pq103_am_both) %>% 
  select(gene_name, protein_name, go, region) %>% kable()
```

| gene\_name | protein\_name                                           | go                                                                   | region |
|:-----------|:--------------------------------------------------------|:---------------------------------------------------------------------|:-------|
| pQBR0006   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0007   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0008   | Putative transmembrane protein                          | <GO:0016021>                                                         | N/A    |
| pQBR0009   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0012   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0014   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0018   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0029   | Putative transmembrane protein                          | <GO:0016021>                                                         | pil    |
| pQBR0030   | Uncharacterized protein                                 |                                                                      | pil    |
| pQBR0032   | Putative pilus assembly-related protein                 |                                                                      | pil    |
| pQBR0033   | Putative pilus biogenesis-like protein                  | <GO:0016021>                                                         | pil    |
| pQBR0034   | Putative transmembrane secretion-related protein        | <GO:0016021>                                                         | pil    |
| pQBR0035   | Putative transmembrane pilus-related protein            | <GO:0016021>                                                         | pil    |
| pQBR0040   | Uncharacterized protein                                 |                                                                      | sam    |
| pQBR0041   | Putative transmembrane protein                          | <GO:0016021>                                                         | sam    |
| pQBR0054   | Uncharacterized protein                                 | <GO:0004527>                                                         | N/A    |
| pQBR0055   | Putative UvrD helicase                                  | <GO:0003677>; <GO:0004003>; <GO:0004527>; <GO:0005524>; <GO:0006281> | N/A    |
| pQBR0059   | Uncharacterized protein                                 | <GO:0035438>                                                         | N/A    |
| pQBR0115   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0150   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0157   | Putative UV resistance protein                          |                                                                      | uvr    |
| pQBR0171   | Ribonuclease HII (RNase HII) (EC 3.1.26.4)              | <GO:0003723>; <GO:0004523>; <GO:0005737>; <GO:0006401>; <GO:0030145> | tra    |
| pQBR0172   | Uncharacterized protein                                 |                                                                      | tra    |
| pQBR0176   | Uncharacterized protein                                 |                                                                      | tra    |
| pQBR0177   | Putative transmembrane protein                          | <GO:0016021>                                                         | tra    |
| pQBR0184   | Uncharacterized protein                                 |                                                                      | tra    |
| pQBR0186   | Putative transmembrane protein                          | <GO:0016021>                                                         | tra    |
| pQBR0187   | Putative exported protein                               |                                                                      | tra    |
| pQBR0189   | Putative conjugative transfer protein                   |                                                                      | tra    |
| pQBR0190   | Uncharacterized protein                                 |                                                                      | tra    |
| pQBR0196   | Putative exported protein                               |                                                                      | N/A    |
| pQBR0197   | Putative exported protein                               |                                                                      | N/A    |
| pQBR0203   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0204   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0210   | Uncharacterized protein                                 | <GO:0016021>                                                         | N/A    |
| pQBR0211   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0227   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0228   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0243   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0266   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0294   | Putative SOS response-associated peptidase (EC 3.4.-.-) | <GO:0008233>                                                         | N/A    |
| pQBR0303   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0392   | Conserved hypothetical exported protein                 |                                                                      | N/A    |
| pQBR0393   | Uncharacterized protein                                 |                                                                      | N/A    |
| pQBR0404   | Uncharacterized protein                                 | <GO:0009279>; <GO:0016021>                                           | N/A    |
| pQBR0447   | Uncharacterized protein                                 |                                                                      | N/A    |

Genes downregulated by ∆gacS specifically for pQBR103 are:

``` r
pq103_gacs_down <- tab_pq103_am %>% filter(sig=="dgacS only") %>% select(gene) %>% pull() %>% unique()

pq103_gene_info %>% filter(gene_name %in% pq103_gacs_down) %>% 
  select(gene_name, protein_name, go, region) %>% kable()
```

| gene\_name | protein\_name                               | go                         | region |
|:-----------|:--------------------------------------------|:---------------------------|:-------|
| pQBR0001   | Putative partitioning protein               |                            | par    |
| pQBR0002   | Putative partitioning protein               | <GO:0003677>               | par    |
| pQBR0017   | Uncharacterized protein                     |                            | N/A    |
| pQBR0019   | Putative transmembrane anchored protein     | <GO:0016021>               | N/A    |
| pQBR0020   | Uncharacterized protein                     |                            | N/A    |
| pQBR0026   | Uncharacterized protein                     |                            | N/A    |
| pQBR0043   | Uncharacterized protein                     | <GO:0003824>               | sam    |
| pQBR0064   | Conserved hypothetical exported protein     |                            | N/A    |
| pQBR0068   | Uncharacterized protein                     |                            | N/A    |
| pQBR0077   | Uncharacterized protein                     |                            | N/A    |
| pQBR0127   | Uncharacterized protein                     |                            | N/A    |
| pQBR0131   | Putative transmembrane protein              | <GO:0016021>               | N/A    |
| pQBR0132   | Uncharacterized protein                     |                            | N/A    |
| pQBR0143   | Uncharacterized protein                     |                            | N/A    |
| pQBR0148   | Uncharacterized protein                     |                            | N/A    |
| pQBR0158   | Putative DNA repair protein ruvB            | <GO:0003684>; <GO:0006281> | uvr    |
| pQBR0179   | Uncharacterized protein                     |                            | tra    |
| pQBR0221   | Conserved hypothetical DNA-binding protein  | <GO:0003677>               | N/A    |
| pQBR0229   | Uncharacterized protein                     |                            | N/A    |
| pQBR0265   | Uncharacterized protein                     |                            | N/A    |
| pQBR0271   | Uncharacterized protein                     |                            | N/A    |
| pQBR0282   | Uncharacterized protein                     |                            | N/A    |
| pQBR0287   | Uncharacterized protein                     |                            | N/A    |
| pQBR0296   | Putative transmembrane protein              | <GO:0016021>               | N/A    |
| pQBR0302   | Putative transmembrane protein              | <GO:0016021>               | N/A    |
| pQBR0307   | Putative acetyltransferase                  | <GO:0016740>               | N/A    |
| pQBR0311   | Putative stringent starvation protein A     |                            | N/A    |
| pQBR0314   | Uncharacterized protein                     |                            | N/A    |
| pQBR0315   | Uncharacterized protein                     |                            | N/A    |
| pQBR0316   | Putative transmembrane protein              | <GO:0016021>               | N/A    |
| pQBR0327   | Putative phage-related protein              |                            | N/A    |
| pQBR0343   | Uncharacterized protein                     |                            | N/A    |
| pQBR0372   | Uncharacterized protein                     |                            | N/A    |
| pQBR0375   | Putative phage-related protein              |                            | N/A    |
| pQBR0412   | Conserved hypothetical HD domain protein    |                            | N/A    |
| pQBR0414   | Uncharacterized protein                     |                            | N/A    |
| pQBR0415   | Uncharacterized protein                     |                            | N/A    |
| pQBR0416   | Conserved hypothetical exported protein     |                            | N/A    |
| pQBR0420   | Uncharacterized protein                     |                            | N/A    |
| pQBR0446   | Uncharacterized protein                     |                            | N/A    |
| pQBR0472   | Putative CheR chemotaxis signalling protein | <GO:0008757>               | che    |

Make a heatmap.

``` r
pq103_de_am_plot <- filter(pq103, gene %in% pq103_de_am_list) %>% 
  left_join(select(tab_pq103_am, region, colour, gene, sig), by="gene")

pq103_de_am_plot_gene_names <- data.frame(pq103_de_am_plot[,c("gene","sig","region")],
                                         logFC=0,
                                         amelioration="region")

pq103_de_am_plot_gene_names$region <- factor(pq103_de_am_plot_gene_names$region,
                                            levels=regions_cols$region)

(plot_fig15 <- ggplot(data=pq103_de_am_plot,
       aes(y=amelioration, x=fct_reorder(gene, -logFC))) +
  geom_tile(colour="black", aes(fill=logFC, size=ifelse(FDR<0.05,"sig","ns"))) + 
  scale_size_manual(values=c(0,0.1), guide=FALSE) +
  scale_fill_gradient2(high="red",mid="grey90",low="blue",limits=c(-3,3)) +
  geom_point(data=pq103_de_am_plot_gene_names, size=0.8,
             aes(colour=region, alpha=region, shape=region)) +
  scale_alpha_manual(values=c(rep(1,6),0), name="") +
  scale_colour_manual(values=regions_cols$colour, name="") +
  scale_shape_manual(values=c(rep(20,6),32), name="") +
  facet_grid(amelioration~sig, scales="free", space="free") + 
        scale_y_discrete(breaks=c("dPFLU4242","PQBR57_0059_V100A","dgacS"),
                     labels=c("∆PFLU4242","PQBR57_0059_V100A","∆gacS")) +
  theme_pub() +
  labs(x="", y="") +
  scale_x_discrete(expand = expansion(add=1)) +
  theme(legend.position="bottom", strip.text=element_blank(),
        panel.border=element_blank(), axis.line=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing.x=unit(5, units="pt")))
```

![](COMPMUT_RNAseq_6c_analysis_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

Load up the pQBR57 image and output both together.

``` r
library(patchwork)

plot_fig14b <- readRDS("../plots/Fig14b.rds")
plot_fig14b <- plot_fig14b + theme_pub() + 
  theme(strip.text=element_blank(),
        panel.border=element_blank(), axis.line=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing.x=unit(5, units="pt"))
plot_fig15 <- plot_fig15 + theme_pub() + 
  theme(legend.position="bottom", strip.text=element_blank(),
        panel.border=element_blank(), axis.line=element_blank(),
        axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        panel.spacing.x=unit(5, units="pt"))

svglite::svglite(height=4, width=3.6, file = "../plots/Fig14_15.svg")
plot_fig14b / plot_fig15 + plot_layout(heights=c(1.4,1))
dev.off()
```

    ## quartz_off_screen 
    ##                 2

**Summary: pQBR103 expression with amelioration follows a very different
pattern to pQBR103, with more genes in general being upregulated than
downregulated. Note that this analysis was using mapping to chromosomal
genes as a measure of correcting for library size, so this may represent
a general underexpression in the wild-type.**

**Removing PFLU4242 seems to upregulate more genes on pQBR103 than
knocking out gacS.**

------------------------------------------------------------------------

**[Back to index.](COMPMUT_index.md)**
