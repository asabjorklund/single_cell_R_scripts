
# creates nD artifical doublets from a collection of cells
# with number of cells in each sample in nCells
# returns the mean number of within sample doublets
# and across sample doublets after nIter randomizations.

simulate_doublets = function(nD, nCells, nIter=100){
  predCells = rep(names(nCells),times=nCells)
  out = c()
  for ( i in 1:nIter ){
    randCells = predCells[sample(1:length(predCells))]
    # pair every 2 cells into 1, list as across,or within sample doublet
    doublets = randCells[1:(nD*2)]
    m = matrix(doublets, ncol = 2)
    same = m[,1] == m[,2]
    hashedD = sum(!same)
    t = table(factor(m[same,1], levels = names(nCells)))
    t["HashDoublet"] = hashedD
    out = rbind(out,t)
  }
  return(colMeans(out))
}


# runs through a range of doublet rates and runs simulation of doublets
# nCells - a named vector with number of hashed cells per hashtag.
# nMulti - total number of double hashtags observed.

# Identifies the doublet rate closest to the observed nMulti and returns
# a data frame with number of doublets in each sample and double-tagged.

predict_hashed_doublet_rates = function(nCells, nMulti, show.plot = TRUE) { 
  range = nMulti:(nMulti*2)
  if (is.null(names(nCells))){
    names(nCells) = LETTERS[1:length(nCells)]
  }
  
  # simulate doublets
  pred = sapply(range, simulate_doublets, nCells)
  
  # find point closest to true nMulti and draw a line there
  pred = data.frame(t(pred))
  pred$totalDoublets = range
  cut.idx = which.min(abs(pred$HashDoublet - nMulti))
  cut = pred$totalDoublets[cut.idx]
  
  # plot results
  if (show.plot){
    p = ggplot(melt(pred, id.vars = "totalDoublets"), aes(x=totalDoublets,y=value, col=variable)) + 
      geom_point() + geom_vline(xintercept = cut) + ggtitle("Simulated doublets")
    print(p)
  }
  
  out = data.frame(nDoublets = t(pred[cut.idx,]), Percentage = t(pred[cut.idx,]/sum(nCells)*100))
  colnames(out) = c("nDoublets","PercentageTotal")
  out$PercentageSample = c(t(pred[cut.idx,1:length(nCells)]/nCells*100), NA,NA)
  
  return(out)
}