################################################################################
### mobileRNA project: sRNA simulated data - mobileRNA method, step 2  ###
################################################################################

library("GenomicRanges")
library("utils")
library("Repitools")

####### TomTom_vs_TomEgg  #########
# location of shortstack folders
path_to_aligment_files <- file.path("../3_de_novo_detection/")

# sample names matching folders
samples <- c("TomTom_1", "TomTom_2", "TomTom_3","TomEgg_1", "TomEgg_2", "TomEgg_3", "TomEgg_4")


gff_alignment <- GenomicRanges::GRangesList()
for (i in samples) {
  file_path <- file.path(path_to_aligment_files, i, "Results.gff3")
  if (file.exists(file_path)) {
    gff_alignment[[i]] <- rtracklayer::import.gff(file_path)
  } else{
    stop("File does not exist:", file_path, "\n")
  }
}
gff_merged <- GenomicRanges::reduce(unlist(gff_alignment),
                                    ignore.strand = TRUE)
gff_merged <- Repitools::annoGR2DF(gff_merged)
locifile_txt <- data.frame(Locus = paste0(gff_merged$chr, ":",
                                          gff_merged$start,"-",
                                          gff_merged$end),
                           Cluster = paste0("cluster_",
                                            seq_len(nrow(gff_merged))))

loci_out <- file.path("../4_sRNA_loci/TomTom_vs_TomEgg_locifile.txt")
utils::write.table(locifile_txt, file = loci_out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)


####### EggEgg_vs_EggTom  #########
# location of shortstack folders
path_to_aligment_files <- file.path("../3_de_novo_detection/")

# sample names matching folders
samples <- c("EggEgg_1", "EggEgg_2", "EggEgg_3","EggEgg_4", "EggTom_1", "EggTom_2", "EggTom_3")


gff_alignment <- GenomicRanges::GRangesList()
for (i in samples) {
  file_path <- file.path(path_to_aligment_files, i, "Results.gff3")
  if (file.exists(file_path)) {
    gff_alignment[[i]] <- rtracklayer::import.gff(file_path)
  } else{
    stop("File does not exist:", file_path, "\n")
  }
}
gff_merged <- GenomicRanges::reduce(unlist(gff_alignment),
                                    ignore.strand = TRUE)
gff_merged <- Repitools::annoGR2DF(gff_merged)
locifile_txt <- data.frame(Locus = paste0(gff_merged$chr, ":",
                                          gff_merged$start,"-",
                                          gff_merged$end),
                           Cluster = paste0("cluster_",
                                            seq_len(nrow(gff_merged))))

loci_out <- file.path("../4_sRNA_loci/EggEgg_vs_EggTom_locifile.txt")
utils::write.table(locifile_txt, file = loci_out, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
