generatePermutation <- function (noTraits, traitIds, numFolds){
    
    numObs = length(noTraits);
    permutMatV = matrix(numeric(0), nrow = numObs, ncol = numFolds);
    permutMatTest = matrix(numeric(0), nrow = numObs, ncol = numFolds);
    
    #generate all test, and validation
    for (id in 1 : numObs) {
        permt <- c()
        permv <- c()
        
        while (length(permt) < numFolds) {
            traitId = traitIds[[id]];
            len = length(traitId);
            if(len == 1) {
                permt = rep(traitId, numFolds);
                permv = rep(traitId, numFolds);
                break
            }	else {
                ranPermT = sample(traitId, len, replace=FALSE);
                ranPermV = c(ranPermT[2: len], ranPermT[1]);
                
                permt = c(permt, ranPermT);
                permv = c(permv, ranPermV );
            }
            
        }
        
       	permutMatTest[id,] = permt[1:numFolds];
        permutMatV[id,] = permv[1:numFolds];
    }
    return(list(permutMatTest, permutMatV));
    
}




