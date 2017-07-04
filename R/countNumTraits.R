countNumTraits <- function(X, flagNan){
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #find number of traits for each observation
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numRows = nrow(X);
    numCols = ncol(X);
    noTraits = rep(0, times = numRows);
    
    traitIds = vector('list', numRows);
    
    for (plantId in 1 : numRows){
        if(flagNan){
            traitId = which(!is.na(X[plantId, ]));
        }
        else{
            traitId = which(X[plantId, ] != 0);
        }
        
        traitIds[[plantId]] = traitId;
        noTraits[plantId] = length(traitId);
    }
    
    indObs = which(noTraits > 0);
    
    return(list(noTraits ,indObs ,traitIds))
    
}
