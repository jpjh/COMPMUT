COMPMUT Sequences 1: Presentation of previously-identified compensatory
mutations
================
jpjh
compiled Mar 2021, edited Jul 2021

[Now published in PLoS
Biology](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001225):

Hall, J. P. J., Wright, R. C. T., Harrison, E., Muddiman, K. J., Jamie
Wood, A., Paterson, S., & Brockhurst, M. A. (2021). Plasmid fitness
costs are caused by specific genetic conflicts enabling resolution by
compensatory mutation. *PLoS Biology*, *19*(10), e3001225.
<https://doi.org/10.1371/journal.pbio.3001225>

**[Back to index.](COMPMUT_index.md)**

------------------------------------------------------------------------

## PFLU4242

Get mutations in PFLU4242 from the following studies:

-   [Carrilero et
    al. 2021](https://journals.asm.org/doi/full/10.1128/mBio.00558-21):
    *P. fluorescens* with pQBR103 and/or pQBR57.
-   [Hall et
    al. 2019](https://www.microbiologyresearch.org/content/journal/micro/10.1099/mic.0.000862):
    *P. fluroescens* with pQBR55.
-   [Hall et
    al. 2017](https://www.nature.com/articles/s41559-017-0250-3): *P.
    fluorescens* and *P. putida* with pQBR57.
-   [Harrison et
    al. 2017](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14080):
    *P. fluorescens* with pQBR103 and bacteriophage phi-2.
-   [Harrison et
    al. 2015](https://linkinghub.elsevier.com/retrieve/pii/S0960982215007186):
    *P. fluorescens* with pQBR103.

#### Carrilero et al. 2021

``` r
comp_muts <- c("PFLU_3777", "PFLU_2189", "PFLU_4242")

carrilero_variants <- read.csv("../ext/carrilero_all.variants.csv")
carrilero_compmuts <- carrilero_variants %>% filter(GeneID %in% comp_muts) %>%
  mutate(position_na = ifelse(GeneID == "PFLU_3777",
                              1 + POS - 4172520,
                              ifelse(GeneID == "PFLU_4242",
                                            1 - (POS - 4685929),
                                            ifelse(GeneID == "PFLU_2189",
                                                   1 + POS - 2372682, "ERROR"))))
```

Output table for manual conversion of amino acids into single letter
code.

``` r
write.table(select(carrilero_compmuts, Clone, GeneName, Change, REF, ALT, position_na),
            file="../ext/carrilero_compmuts.csv", sep=",",
            quote=FALSE, row.names=FALSE)
```

Edit manually to get single-letter codes.

Reduce to get unique mutations per population.

``` r
carr_compmuts <- read.csv("../ext/carrilero_compmuts_edit.csv") %>% 
  mutate(pop = gsub("-.*", "", Clone),
         position_aa = as.integer(gsub("(.*)_.*", "\\1",
                         gsub("[^0-9_]", "", SingleLetter))),
         study = "Carrilero et al. (2021)",
         target = ifelse(GeneName == "PFLU_3777", "gacS",
                         ifelse(GeneName == "PFLU_4242", "PFLU4242", 
                                ifelse(GeneName=="gacA", "gacA", "error"))),
         trt = str_sub(Clone, 4, 4),
         plasmid = ifelse(trt %in% c(1,5), "pQBR103",
                          ifelse(trt %in% c(2,6), "pQBR57",
                                 ifelse(trt %in% c(4,8), "pQBR103 and pQBR57", "ERROR")))) %>%
  rename(clone = pop, mutation_aa = SingleLetter) %>%
  select(study, clone, target, mutation_aa, position_aa, position_na, plasmid) %>% unique()
```

#### Hall et al. 2019

``` r
pq55_variants <- read.csv("../ext/TabS6_GeneticsData.csv")

pq55_compmuts <- pq55_variants %>% 
  pivot_longer(cols = c(gacS, gacA, PFLU4242, rpo), names_to = "target", values_to = "mutation") %>%
  filter(!(mutation %in% c("ND","wt")) & summary != "rpo") %>%
  filter(target != "del") %>%
  mutate(position_aa = as.integer(gsub("[^0-9_]", "", mutation)),
         position_na = as.integer(gsub("[^0-9_]", "", mutation_na)),
         study = "Hall et al. (2019)",
         lineage = as.character(lineage), plasmid="pQBR55") %>%
  rename(clone = lineage, mutation_aa = mutation) %>% select(!summary)
```

#### Hall et al. 2017

``` r
soili_compmuts <- read.csv("../ext/compmuts_edit.csv") %>%
  filter(target != "par 5'") %>%
  mutate(clone = paste(sample_number, cul, mer, spc, rep, sep="-"),
         position_na = mutation_pos,
         position_aa = as.integer(gsub("[^0-9_]", "", mutation_aa)),
         mutation_aa = ifelse(mutation_aa == "0",
                              paste(type,"ins",sep=""), mutation_aa),
         study = "Hall et al. (2017)", plasmid="pQBR57") %>%
  select(study, clone, target, position_aa, position_na, mutation_aa, plasmid)
```

#### Harrison et al. 2017

``` r
plapha_compmuts <- read.csv("../ext/Harrison2017_data.csv") %>%
  filter(LOCUS_ID %in% c("PFLU_4242","PFLU_3777","PFLU_2189")) %>%
  mutate(position_na = ifelse(LOCUS_ID == "PFLU_3777",
                                     1 + POSITION - 4172520,
                                     ifelse(LOCUS_ID=="PFLU_4242",
                                            1 - (POSITION - 4685929),
                                            ifelse(LOCUS_ID=="PFLU_2189",
                                                   1 + POSITION - 2372682, "ERROR"))),
         target = factor(GENE_NAME..IF.AVALIBABLE., 
                         levels=c("gacA","gacS","PFLU_4242"),
                         labels=c("gacA","gacS","PFLU4242")),
         plasmid = "pQBR103") %>%
  rename(ref = WILD.TYPE.BASE.SEQUENCE, alt = MUTATION.BASE.SEQUENCE, 
         clone = SAMPLE.ID..Library.code...population.ID...extraction.number.) %>%
  select(clone, ref, alt, target, position_na, plasmid)
```

Need to annotate these!

``` r
write.table(plapha_compmuts,
            file="../ext/plapha_compmuts.csv", sep=",",
            quote=FALSE, row.names=FALSE)
```

Annotated manually.

``` r
pp_compmuts <- read.csv("../ext/plapha_compmuts_edits.csv") %>%
  mutate(study = "Harrison et al. (2017)") %>% select(!(c(alt, ref)))
```

#### Harrison et al. 2015

``` r
parr_compmuts <- read.csv("../ext/Re-seq master mutations list.csv") %>%
  filter(ID %in% c("PFLU4242", "PFLU3777", "gacA")) %>%
    mutate(position_na = as.integer(ifelse(ID == "PFLU3777",
                                     1 + position - 4172520,
                                     ifelse(ID=="PFLU4242",
                                            1 - (position - 4685929),
                                            ifelse(ID=="gacA",
                                                   1 + position - 2372682, "ERROR")))),
           mutation_aa = sapply(strsplit(INFO, "\\|"), `[`, 3),
           target = ifelse(ID == "PFLU3777",
                           "gacS", ifelse(ID=="PFLU4242",
                                  "PFLU4242", ifelse(ID=="gacA",
                                         "gacA", "ERROR"))),
           mutation_na = paste(Ref, position_na, Difference, sep=""),
           plasmid = "pQBR103",
           study = "Harrison et al. (2015)") %>%
  rename(clone = Clone) %>%
  select(clone, target, position_na, mutation_aa, mutation_na, target, plasmid, study)
```

One entry, P40-5, has a deletion from Leu437 to Gln442, and the amino
acid mutations haven’t copied across for P32-5. Replace this manually.

``` r
parr_compmuts[25,"mutation_aa"] <- "L437_Q442del"
parr_compmuts[20,"mutation_aa"] <- "V148M"
```

### Bind all together and plot.

``` r
compmuts <- bind_rows(parr_compmuts, pp_compmuts, carr_compmuts, pq55_compmuts, soili_compmuts) %>%
  select(study, plasmid, clone, target, position_na, mutation_na, position_aa, mutation_aa)

write.csv(compmuts, file = "../data/COMPMUT_mutations.csv", row.names=FALSE)
```

``` r
ggplot(data=compmuts, aes(x=position_na, y=plasmid, colour=study)) + geom_point() + 
  geom_text(aes(label=mutation_aa), angle=90, hjust=0) +
  facet_grid(. ~ target, scales = "free_x", space = "free_x")
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 rows containing missing values (geom_text).

![](COMPMUT_seqs_1_compmuts_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

OK. Needs to be tidied up, but this is the basics.

Load up annotation of domains.

``` r
domains <- read.csv("../ext/target_domains.csv") %>% 
  mutate(lab_pos = ((end_pos - start_pos)/2) + start_pos)
```

Tidy up annotations to just have summary. Interested readers can consult
the supplementary information.

``` r
genelengths <- data.frame(target = c("PFLU4242","PQBR57_0059","gacS","gacA"),
                          startpos = c(1,-25,1,1),
                          endpos = c(1584,669,2754,642))

compmuts$mut_type <- "missense"
compmuts[grep("fs", compmuts$mutation_aa),"mut_type"] <- "frameshift"
compmuts[grep("del|dup", compmuts$mutation_aa),"mut_type"] <- "deletion/duplication"
compmuts[grep("ins", compmuts$mutation_aa),"mut_type"] <- "insertion"
compmuts[grep("\\*", compmuts$mutation_aa),"mut_type"] <- "stop"
```

Summarise the table to get counts for each mutation.

``` r
compmuts_summ <- compmuts %>% 
  group_by(target, position_na, mutation_aa, plasmid, mut_type) %>%
  summarise(count = n(),
            counts = ifelse(count>1, ">1", "1"))
```

    ## `summarise()` has grouped output by 'target', 'position_na', 'mutation_aa', 'plasmid'. You can override using the `.groups` argument.

Plot.

``` r
pd <- position_dodge(width=0.3)

plot_mutations <- function(t){
  muts <- filter(compmuts_summ, target==t)
  annot <- filter(domains, target==t)
  lengths <- filter(genelengths, target==t)
  x_scale <- c(1,seq(500,lengths$endpos,500),lengths$endpos)
  
  ggplot(data=muts) +
    geom_segment(aes(y=position_na, yend=position_na, x=plasmid, 
                     group=plasmid), 
                 position=pd, xend=0, size=0.2) +
    geom_point(aes(y=position_na, x=plasmid, shape=mut_type, group=plasmid, colour=plasmid, size=counts),
               position=pd) +
    geom_rect(data=annot,
              aes(ymax = end_pos, ymin = start_pos, fill = fill_col),
              size = 0, xmin = 0, xmax = 0.5) +
    geom_text(data=annot, size=2,
              aes(y = lab_pos, label = domain_name), 
              colour = "white", vjust = 0.5, hjust = 0.5, x = 0.25) +
    scale_y_continuous(limits=c(lengths$startpos, lengths$endpos),
                       name="",
                       breaks=x_scale) +
    scale_size_manual(breaks=c("1",">1"), values=c(1.5,4.5)) + 
    # use size of 1 when just 1, use size of 2 for >1 parallel mutation
    scale_shape_manual(breaks=c("missense","frameshift","deletion/duplication","insertion","stop"),
                       values=c("\u25CF", "\u25A0", "\u25B2", "\u25BC", "\u2736"), name = "") +
    scale_colour_manual(breaks=c("pQBR55","pQBR57","pQBR103","pQBR103 and pQBR57"),
                        values=c("hotpink","springgreen3","dodgerblue","goldenrod2"), 
                        name = "") +
    scale_fill_manual(values = annot$fill_col, guide = "none") +
    coord_flip(xlim=c(0,5), expand=FALSE) +
    theme_minimal(base_size = 6) +
    theme(panel.grid=element_blank(), axis.text.y = element_blank(),
          axis.line.x=element_line(size=0.4), axis.ticks.x=element_line(size=0.4), 
          axis.text.x=element_text(colour="black"), axis.title.y = element_blank(),
          legend.position = "none")
}
```

Output as .svg.

Output to scale by dividing the gene lengths by a set factor, so the
scale images can be compiled in Inkscape.

``` r
scalefactor <- 800

output_svg <- function(x) {
  gene <- x
  w <- genelengths[genelengths$target==gene,"endpos"]/scalefactor
  filename <- paste("../plots/Fig_mutations/Fig_mutations_",gene,".svg", sep="")
  svglite::svglite(height=1, 
                 width=w, 
                 file=filename)
  print(plot_mutations(gene))
  dev.off()
}

output_svg("PFLU4242")
```

    ## Warning: Removed 1 rows containing missing values (geom_segment).

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## quartz_off_screen 
    ##                 2

``` r
output_svg("PQBR57_0059")
```

    ## quartz_off_screen 
    ##                 2

``` r
output_svg("gacS")
```

    ## quartz_off_screen 
    ##                 2

``` r
output_svg("gacA")
```

    ## quartz_off_screen 
    ##                 2

Create a legend.

``` r
svglite::svglite(height=10, 
                 width=7, 
                 file="../plots/Fig_mutations/Legend.svg")
plot_mutations("gacS") + theme(legend.position = "right")
dev.off()
```

    ## quartz_off_screen 
    ##                 2

Put all together in Inkscape.

------------------------------------------------------------------------

**[Back to index.](COMPMUT_index.md)**
