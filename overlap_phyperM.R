overlap_phyperM<-function(v1,v2,bg=length(v1), plot=FALSE,
                          title="Overlap", remove.diag = FALSE, with.total = TRUE,
                          to.file = NULL, silent=FALSE, nsize = 15){
   # takes two vectors of metadata and creates a matrix with pvalue for overlap between elements
   # phyper test uses all entries in L as background if bg is not specified.
   if (length(v1) != length(v2)){ stop("Error, the vectors are not of the same length")}
   L1 =   split(1:length(v1), v1)
   L2 = split(1:length(v2), v2)
  
   nL1<-length(L1)
   nL2<-length(L2)
   M<-mat.or.vec(nL1,nL2)
   P<-mat.or.vec(nL1,nL2)
   P[,]<-1
   for (i in 1:nL1){
      for (j in  1:nL2){
        M[i,j]<-length(intersect(L1[[i]],L2[[j]]))
        if (M[i,j] ==0) {
          P[i,j] <- 1
        }else{
          P[i,j]<-phyper(M[i,j]-1,length(L1[[i]]),bg-length(L1[[i]]),length(L2[[j]]), lower.tail = FALSE)
        }
      }
   }
   colnames(P)<-names(L2)
   rownames(P)<-names(L1)
   colnames(M)<-names(L2)
   rownames(M)<-names(L1)

   P[M==0]<-1
   # still may have values that are zero
   pseudo = min(P[P>0])*0.1

   if(remove.diag){
      diag(P) = NA
   }

   # ad column/row with total genes per list.
   # lower square is total unique genes in the 2 lists.
   if (with.total){
      nTotal1 = unlist(lapply(L1,length))
      nTotal2 = c(unlist(lapply(L2,length)),length(unique(c(unlist(L1),unlist(L2)))))
      M = cbind(M,nTotal1)
      M = rbind(M,nTotal2)
      P = cbind(P, rep(NA,nrow(P)))
      P = rbind(P, rep(NA,ncol(P)))

   }
   suppressMessages(require(pheatmap))
   if (is.null(to.file)){
      pl = pheatmap(-log10(P+pseudo), cluster_rows = F, cluster_cols = F, display_numbers = M,
                    main = title, silent = silent, fontsize_number = nsize)
   }else{
      pheatmap(-log10(P+pseudo), cluster_rows = F, cluster_cols = F, display_numbers = M,
               main = title, filename = to.file, silent = TRUE, fontsize_number = nsize)
   }
   return(list(P=P,M=M, plot=pl))
}