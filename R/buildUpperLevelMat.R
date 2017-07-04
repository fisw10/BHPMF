buildUpperLevelMat <- function(input_dir, cvInd, traitNames, num_hierarchy_levels, numFolds, verbose) {
  
  testFlag <- 0
  hierarchy.info.idx <- NULL	
  level.names <- list()
  file.name <- paste(input_dir, '/processed_hierarchy_info.Rda', sep = "");
  load(file.name);
  
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #create upper level Matrices
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for (datasetId in 1 : numFolds) {
    
    if(verbose) {
      cat(datasetId)
    }
    if(cvInd > 0) {
      dir_name = paste(input_dir, "/fold", datasetId, "/Tuning/fold", cvInd, sep = "")
      secondCV = cvInd;
    } else {
      dir_name = paste(input_dir, "/fold", datasetId, sep = "")
      secondCV = 0;
    }
    
    file_name <- paste(dir_name, "/Ytrain", num_hierarchy_levels, ".txt", sep = "")
    tr <- as.matrix(read.table(file_name, sep="\t"))
    X <- as.matrix(sparseMatrix(tr[,1], tr[,2], x = tr[,3]));
    rm(tr)
    
    num_cols <- ncol(X);
    for (level in 1 : (num_hierarchy_levels-1)){
      
      msg = paste("building upper level matrix at level ", level, "\n")
      if(verbose) {
        cat(msg)
      }
      #%%%%%%%%%%%%%%%%%
      #create upper level matrix at level
      #%%%%%%%%%%%%%%%%%
      numSpecies <- length(level.names[[num_hierarchy_levels - level]])
      #	    treeinds_eachspecies <- vector('list', numSpecies)
      
      XSpecies <- matrix(data=NA, numSpecies, num_cols)
      for (ind in 1 : numSpecies) {
        
        treeinds_eachspecies <- which(hierarchy.info.idx[, level] == ind)
        
        info <- X[treeinds_eachspecies, ];
        
        if(length(treeinds_eachspecies) == 1) {
          info = t(as.matrix(info, nrow = 1, ncol = length(info)))
        }
        
        for (trait in 1: num_cols){
          ss <- sum(info[, trait])
          len <- length(which(info[, trait] != 0))
          XSpecies[ind, trait] <- ss / len   
        }    
      }
      
      # XSpecies[is.nan(XSpecies)] <- 0
      
      #%%%%%%%%%%%%%%%%%%%%%%%%
      # write datasets
      #%%%%%%%%%%%%%%%%%%%%%%%%
      
      find_folds(XSpecies, num_hierarchy_levels-level, 1, -1, dir_name)
      
      # #find non-missing elements
      # Idx <- which(!is.na(XSpecies), arr.ind=TRUE);
      # val <- XSpecies[Idx];
      
      # #Y: 3-column matrix: sparse representation of trait info
      # #each row is a triple (i, j, value)
      # #i: observation index (row index)
      # #j: trait index (column index)
      # #vlaue: the observed value
      # Y <- cbind(Idx, val);
      # rm(Idx, val);
      
      # file_name <- paste(dir_name, "/Ytrain", num_hierarchy_levels-level, ".txt", sep = "")
      # write.table(Y, file=file_name, sep="\t", col.names = F, row.names = F)
    }
  }
}