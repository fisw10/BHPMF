create_upper_level_matrix <- function (X, phyinfo, level, numSpecies) {

	#%%%%%%%%%%%%%%%%%
    #create upper level matrix at level
    #%%%%%%%%%%%%%%%%%
	num_cols <- ncol(X)
    X_upp_level <- matrix(numeric(0), numSpecies, num_cols);
    for (ind in 1 : numSpecies) {
	    treeinds_eachspecies <- which(phyinfo[, level] == ind);

        info <- X[treeinds_eachspecies, ];

        if(length(treeinds_eachspecies) == 1)
    	    info = t(as.matrix(info, nrow = 1, ncol = length(info)));

		for (trait in 1: num_cols) {
        	ss <- sum(info[,trait]);
            len <- length(which(info[,trait] != 0));
            X_upp_level[ind, trait] <- ss / len;
        }
	}
    X_upp_level[is.nan(X_upp_level)] <- 0;
	return(X_upp_level);
}