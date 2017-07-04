find_folds <- function(X, level, num_folds, secondCV, output_dir) {
    res <- countNumTraits(X, 0);
    num_available_features_per_row <- res[[1]];
    indObs <- res[[2]];
    traitIds <- res[[3]];
    rm(res)
    
    res <- generatePermutation(num_available_features_per_row, traitIds, num_folds);
    permutMatTest <- res[[1]];
    rm(res)
    
    # find index of rows with more than two avialble features
    ii2 = which(num_available_features_per_row >= 2);
    jj_test = matrix(permutMatTest[ii2, ], ncol = ncol(permutMatTest));
    ii_test = matrix(rep(ii2, num_folds), nrow = length(ii2), ncol = num_folds);
    
    if(secondCV >= 0) {
        for (fold in 1 : num_folds) {
            dir_name <- paste(output_dir, "/fold", fold, sep="")
            # If directry does not exists, create the directory
            if (!file.exists(dir_name)) {
                dir.create(dir_name)
            }
            split_data(X, ii_test[, fold], jj_test[, fold], level, dir_name)
        }
    } else {
        split_data(X, ii_test[, 1], jj_test[, 1], level, output_dir)
    }
}
