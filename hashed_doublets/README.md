# Predict intra-sample doublet rates based on hashing

Given the number of cells per hashtag in an experiment and the number of double tags observed,  the intra-sample doublet rate can be calculated for each sample.

Assuming that the sample distribution in the hashed cells and in the doublets is the same, e.g. that there is no sample bias in doublet creation. Also, assumes that all multiplets consists of 2 cells.

The function in `predict_hashed_doublet_rates.R` will create artificial doublets from the observed distribution, randomly selected 100 times. Simulations are done across different doublet rates and then identifies the rate closest to the observed for across hashtag doublets.

Examples are provided in the [Notebook](doublet_rate_predictor.md)
