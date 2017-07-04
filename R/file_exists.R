CheckPreprocessFilesExist <- function(tmp.dir, num.folds, num.hierarchy.levels, prediction.level) {
	# Check if hierarchy related files exist
	if (!file.exists(paste(tmp.dir, "/num_parents.txt", sep=""))) return(FALSE)
	if (!file.exists(paste(tmp.dir, "/num_children.txt", sep=""))) return(FALSE)
	if (!file.exists(paste(tmp.dir, "/hierarchy_info.txt", sep=""))) return(FALSE)
	if (!file.exists(paste(tmp.dir, "/processed_hierarchy_info.Rda", sep=""))) return(FALSE)

	# Check training files	
	for (fold in 1:num.folds) {
		for (level in 1:num.hierarchy.levels) {
			if (!file.exists(paste(tmp.dir, "/fold", fold, "/Ytrain", level , ".txt", sep=""))) return(FALSE)
		}
	}

	# Check test file
    for (fold in 1:num.folds) {
	    if (!file.exists(paste(tmp.dir, "/fold", fold, "/Ytest", prediction.level, ".txt", sep=""))) return(FALSE)
    }

	return(TRUE)
}

CheckPreprocessFilesExistTuning <- function(tmp.dir, num.folds, used.num.hierarchy.levels, prediction.level) {
    # Check if hierarchy related files exist
    if (!file.exists(paste(tmp.dir, "/num_parents.txt", sep=""))) return(FALSE)
    if (!file.exists(paste(tmp.dir, "/num_children.txt", sep=""))) return(FALSE)
    if (!file.exists(paste(tmp.dir, "/hierarchy_info.txt", sep=""))) return(FALSE)

    # Check training files
    for (fold in 1:num.folds) {
        for (level in 1:used.num.hierarchy.levels) {
            if (!file.exists(paste(tmp.dir, "/fold", fold, "/Ytrain", level , ".txt", sep=""))) return(FALSE)
        }
    }

    # Check test file
    for (fold in 1:num.folds) {
        if (!file.exists(paste(tmp.dir, "/fold", fold, "/Ytest", prediction.level, ".txt", sep=""))) return(FALSE)
    }

	return(TRUE)
}