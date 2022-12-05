# copied from:
#https://ucdavis-bioinformatics-training.github.io/2020-Advanced_Single_Cell_RNA_Seq/data_analysis/VDJ_Analysis_fixed
# and modified to parse both heavy and light chain info.

add_clonotype_hl <- function(tcr_prefix, seurat_obj, type="t"){
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))

  
  # remove barcodes that are not in seurat_obj
  tcr = tcr[tcr$barcode %in% colnames(seurat_obj),]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # select columns
  cn =  colnames(tcr)
  cn.chains = cn[c(3,5,6,7,8,9,10,13:28,30)]
  cn.both = cn[c(1,29,31)]
  
  # get top H/L chain per barcode
  tcr$HL = ifelse(tcr$chain == "IGH","H","L")
  cn.chains = c(cn.chains, "n")
  tcr = split(tcr, tcr$barcode)
  
  # if multiple H/L chains, select the one with highest number of UMIs.
  # also flag ones with multiple sequences.
  # for IGL/IGK select the one with highest umi.
  #tcr = lapply(tcr, function(x) x %>% group_by(HL) %>%  add_tally() %>% filter(umis == max(umis)))
  
  tcr = lapply(tcr, function(x) x %>% group_by(HL) %>%  add_tally() %>% slice(which.min(umis)))
  
  # merge all into one row per barcode
  tcr = lapply(tcr, function(x) {
    x = data.frame(x)
    new = x[1,cn.both]
    light = x[x$HL == "L",cn.chains]
    colnames(light) = paste("LC",colnames(light), sep = "_")
    if (nrow(light) == 0){ light[1,] = NA }
    heavy = x[x$HL == "H",cn.chains]
    colnames(heavy) = paste("HC",colnames(heavy), sep = "_")  
    if (nrow(heavy) == 0){ heavy[1,] = NA }
    return(data.frame(cbind(new,light,heavy)))
  })
  
  
  all = Reduce(rbind,tcr)
  rownames(all) = names(tcr)
  
  print("Doublet LC vs HC")
  print(table(all$LC_n, all$HC_n))
  
  # Add to the Seurat object's metadata.
  colnames(all) <- paste(type, colnames(all), sep="_")
  seurat_obj <- AddMetaData(object=seurat_obj, metadata=all)
  return(seurat_obj)
}



add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
    tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))

    # Remove the -1 at the end of each barcode.
  
    tcr <- tcr[!duplicated(tcr$barcode), ]

    # Only keep the barcode and clonotype columns. 
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"

    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))

    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"


    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3)]
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    colnames(tcr) <- paste(type, colnames(tcr), sep="_")
    # Add to the Seurat object's metadata.
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(clono_seurat)
}
