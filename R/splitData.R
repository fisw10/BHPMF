split_data <- function(X, ii_test, jj_test, level, output_dir) {

	num_elem <- nrow(X) * ncol(X)
	subIdx <- matrix(1:num_elem, nrow(X), ncol(X))
	
	# Create the X and Y matrices for test data (Xtest, Yttest)
	indTest <- subIdx[cbind(ii_test, jj_test)]
	Xtest <- matrix(0, nrow(X), ncol(X))
	Xtest[indTest] <- X[indTest]
	Idx <- which(Xtest != 0, arr.ind=TRUE)
	Ytest <- cbind(Idx, Xtest[Idx])
	
	# Create the X and Y matrices for training data (Xtrain, Ytrain)
	indTrain <- setdiff(1:num_elem, indTest)
	Xtrain <- matrix(0, nrow(X), ncol(X))
	Xtrain[indTrain] <- X[indTrain]
	Idx <- which(Xtrain != 0, arr.ind=TRUE)
	Ytrain <- cbind(Idx, Xtrain[Idx])

	# Write training data
	file_name = paste(output_dir, "/Ytrain", level, ".txt", sep = "")
	write.table(Ytrain, file=file_name, sep="\t", col.names = F, row.names = F)

	# Write test data
	file_name = paste(output_dir, "/Ytest", level, ".txt", sep = "")
	write.table(Ytest, file=file_name, sep="\t", col.names = F, row.names = F)
}
