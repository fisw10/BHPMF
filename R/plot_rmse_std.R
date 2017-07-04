PlotRmseVsStd <- function(test.res, test.std, output.plot.path) {
# This function plot RMSE vs Std
  #
  # Args:
  #   test.res: A vector containing the residual of BHPMF.
  #   test.std: A vector containing the std of BHPMF. It is the same size as test.res.
  #   output.plot.path: If present, the output plot will be saved in this path.
  #

	if (length(test.res) != length(test.std)) {
	   stop("test.res should be of the same length as test.std.")
	}

	sort.res <- sort(test.std, index.return = TRUE)
	sorted.std <- sort.res$x
	sorted.res <- test.res[sort.res$ix]
	
	numObs <- length(test.std)
	numBatches <- 10

	sizeBatch <- floor(numObs / numBatches)

	rmse <- rep(0, numBatches)

	for (batch in 1 : numBatches) {
    	beg <- (batch-1)*sizeBatch+1
    	end <- batch*sizeBatch
    	Idx <- beg:end
    	rmse[batch] <- sqrt(mean(sorted.res[Idx]^2))
	}

	if (!missing(output.plot.path)) {
		pdf(output.plot.path)
		plot(rmse, main="", xlab="Percentage of Data with Ascending Std", ylab="Average RMSE", xlim=c(1, 10), xaxt="n", type = "b")
		axis(1, at=seq(1, 10, by = 1), labels = paste(1:10, "0%", sep=""), las=2)
		dev.off()
	} else {
        plot(rmse, main="", xlab="Percentage of Data with Ascending Std", ylab="Average RMSE", xlim=c(1, 10), xaxt="n", type = "b")
        axis(1, at = seq(1, 10, by = 1), labels = paste(1:10, "0%", sep=""), las=2)
	}
}