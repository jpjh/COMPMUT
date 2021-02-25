## Function for making RNAseq data into long format
# Requires reshape2

makeLong <- function(x) {
  tab_long <- reshape2::melt(x[,grep("(.logCPM)|(logFC)|(FDR)|(gene)",colnames(x))], 
                            id.vars=c("gene"))
  tab_long$contrast <- gsub("(.logCPM)|(.logFC)|(.FDR)", "", tab_long$variable)
  tab_long$type <- gsub("^.*\\.", "", tab_long$variable)
  tab_long2 <- reshape2::dcast(tab_long, gene+contrast~type)
  tab_long2$plasmid <- factor(sapply(strsplit(tab_long2$contrast, "\\."), `[`, 1),
                                         levels=c("none","pQBR57","pQBR103"))
  tab_long2$amelioration <- factor(sapply(strsplit(tab_long2$contrast, "\\."), `[`, 2),
                                              levels=c("wt","dPFLU4242","dgacS","PQBR57_0059_V100A"))
  return(tab_long2)
}
