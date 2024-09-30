# visualize overlap of MIF and MIF-AS1

library(Gviz)
library(GenomicRegionTrack)
library(ggplot2)

# chromosome
itrack <- IdeogramTrack(genome = "hg38", chromosome = "22")
itrack@dp@pars[["fontsize"]] <- 24

# genome axis
gtrack <- GenomeAxisTrack()
gtrack@dp@pars[["fontsize"]] <- 20
gtrack@dp@pars[["size"]] <- 1.2

# get gene annotation
chr22 <- data.frame(
  chromosome = rep("chr22",6),
  start = c(23894383, 23894772, 23895040, 23893710, 23895829, 23898844),
  end = c(23894582, 23894944, 23895223, 23895532, 23896110, 23898930),
  width = c((23894582-23894383+1), (23894944-23894772+1), (23895223-23895040+1), (23895532-23893710+1), (23896110-23895829+1), (23898930-23898844+1)),
  strand = c("+", "+", "+", "-", "-", "-"),
  feature = rep("protein_coding",6),
  gene = c("ENSG00000240972", "ENSG00000240972", "ENSG00000240972", "ENSG00000218537", "ENSG00000218537", "ENSG00000218537"),
  exon = rep("NA",6),
  transcript = rep("NA",6),
  symbol = c("MIF   ", "MIF   ", "MIF   ", "MIF-AS1", "MIF-AS1", "MIF-AS1"),
  stringsAsFactors = TRUE
)


# customize gene region tracks
mif <- GeneRegionTrack(head(chr22,3))
mif_as1 <- GeneRegionTrack(tail(chr22,3))

mif@dp@pars[["fill"]] <- "#96C09F"
mif@dp@pars[["fontsize.group"]] <- 24

mif_as1@dp@pars[["fill"]] <- "#548235"
mif_as1@dp@pars[["fontsize.group"]] <- 24

# add empty tracks to control size
empty_1 <- GeneRegionTrack()
empty_1@dp@pars[["size"]] <- 0.3
empty_1@name <- ""
empty_2 <- GeneRegionTrack()
empty_2@dp@pars[["size"]] <- 0.5
empty_2@name <- ""

# plot

png("figures/supplement_4b.png", width=2500, height=800, res=300)

plotTracks(list(itrack, empty_1, gtrack, empty_2, mif_as1, empty_1, mif), from = 23893000, to = 23899000, 
           transcriptAnnotation = "symbol", collapse=FALSE, mergeGroups=FALSE)

dev.off()
