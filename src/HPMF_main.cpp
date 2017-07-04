//
//  HPMF_main.cpp
//
//
//  Created by Farideh Fazayeli on 7/8/13.
//
//

#include <iostream>

#include <R.h>      // R functions
#include <Rmath.h>  // R math
#include <time.h>
#include <Rdefines.h>

#include "HPMF.h"
#include <math.h>
using namespace std;


int XDemo(int numSamples, double *allRMSE, double *allRMSE2, int fold_id,
          double &rmse, SEXP args, const char* mean_file_path,
          const char* std_file_path, SEXP numNodesPerLevelR) {
    
    int num_columns = 13, num_latent_features = 15;
    int num_hierarchy_level = 4, predict_level, gaps, burn;
    int saveFlag = 0;
    int used_num_hierarchy_level;
    int outWholeFlag = 0;
    int verbose;
    
    double *numNodesPerLevel =  REAL(numNodesPerLevelR);
    
    burn = INTEGER(getListElement(args, "Burn"))[0];
    gaps = INTEGER(getListElement(args, "Gaps"))[0];
    saveFlag = INTEGER(getListElement(args, "SaveFileFlag"))[0];
    outWholeFlag = INTEGER(getListElement(args, "OutWholeFlag"))[0];
    num_columns = INTEGER(getListElement(args, "NumTraits"))[0];
    num_latent_features = INTEGER(getListElement(args, "NumFeats"))[0];
    num_hierarchy_level = INTEGER(getListElement(args, "NumHierarchyLevel"))[0];
    used_num_hierarchy_level =
        INTEGER(getListElement(args, "UsedNumHierarchyLevel"))[0];
    predict_level = INTEGER(getListElement(args, "PredictLevel"))[0];
    verbose = INTEGER(getListElement(args, "Verbose"))[0];
    SEXP dirEXP = getListElement(args, "InputDir");
    const char *input_dir = CHAR(STRING_ELT(dirEXP, 0));

    int numPlants;
    GetRNGstate();
    
    HPMF *obj = new HPMF(num_columns, num_latent_features,numNodesPerLevel,
                         predict_level, num_hierarchy_level, used_num_hierarchy_level,
                         gaps, burn, numSamples, fold_id, input_dir, saveFlag, outWholeFlag,
                         mean_file_path, std_file_path, verbose);
    rmse = obj->runBHPMF(allRMSE, numPlants);
    
    PutRNGstate();
    return numPlants;
}

extern "C" {
    SEXP DemoHPMF(SEXP args, SEXP mean_file_pathSexp, SEXP std_file_pathSexp,
                  SEXP numNodesPerLevelR) {
        
        int numSamples = INTEGER(getListElement(args, "NumSamples"))[0];
        int fold_id = INTEGER(getListElement(args, "DatasetId"))[0];
        
        const char *mean_file_path = CHAR(STRING_ELT(mean_file_pathSexp, 0));
        const char *std_file_path = CHAR(STRING_ELT(std_file_pathSexp, 0));
        
        double rmse;
        double *allRMSE = (double *) R_alloc(numSamples, sizeof(double));
        double *allRMSE2 = (double *) R_alloc(numSamples, sizeof(double));
        int numPlants = XDemo(numSamples, allRMSE, allRMSE2, fold_id, rmse, args,
                              mean_file_path, std_file_path, numNodesPerLevelR);
        
        SEXP result, resultNames, rmse_r, numPlants_r, all_rmse_r, all_rmse_r2;
        
        int nProtected = 0;
        PROTECT(all_rmse_r = NEW_NUMERIC(numSamples));
        ++nProtected;
        
        double *allRMSE_r  =  (double *) R_alloc(numSamples, sizeof(double));
        allRMSE_r = NUMERIC_POINTER(all_rmse_r);
        for(int ii = 0; ii < numSamples; ii++) {
            allRMSE_r[ii] = allRMSE[ii];
        }
        
        PROTECT(all_rmse_r2 = NEW_NUMERIC(numSamples));
        ++nProtected;
        
        double *allRMSE_r2  =  (double *) R_alloc(numSamples, sizeof(double));
        allRMSE_r2 = NUMERIC_POINTER(all_rmse_r2);
        for(int ii = 0; ii < numSamples; ii++) {
            allRMSE_r2[ii] = allRMSE2[ii];
        }
        
        PROTECT(numPlants_r = NEW_INTEGER(1));
        ++nProtected;
        int *numP_r = (int *) R_alloc(1, sizeof(int));
        numP_r = INTEGER_POINTER(numPlants_r);
        numP_r[0] = numPlants;
        
        PROTECT(rmse_r = NEW_NUMERIC(1));
        ++nProtected;
        double *RMSE_r = (double *) R_alloc(1, sizeof(double));
        RMSE_r = NUMERIC_POINTER(rmse_r);
        RMSE_r[0] = rmse;
        
        const char *names[4] = {"allRMSE", "numPlants", "RMSE", "allRMSE2"};
        
        PROTECT(resultNames = allocVector(STRSXP, 4));
        ++nProtected;
        
        for(int i = 0; i < 4; i++)
            SET_STRING_ELT(resultNames,i,mkChar(names[i]));
        
        
        PROTECT(result = allocVector(VECSXP, 4));
        ++nProtected;
        
        SET_VECTOR_ELT(result, 0, all_rmse_r);
        SET_VECTOR_ELT(result, 1, numPlants_r);
        SET_VECTOR_ELT(result, 2, rmse_r);
        SET_VECTOR_ELT(result, 3, all_rmse_r2);
        setAttrib(result, R_NamesSymbol, resultNames);
        UNPROTECT(nProtected);
        
        return result;
        
    }
    
    void writeFilledGap(SEXP args, SEXP inName, SEXP outName){
        
        int numSamples = INTEGER(getListElement(args, "numSamples"))[0];
        int num_columns = INTEGER(getListElement(args, "numTraits"))[0];
        int numPlants = INTEGER(getListElement(args, "numPlants"))[0];
        int gap = INTEGER(getListElement(args, "gap"))[0];
        int burn = INTEGER(getListElement(args, "burn"))[0];
        const char *inFilePath = CHAR(STRING_ELT(inName, 0));
        const char *outFilePath = CHAR(STRING_ELT(outName, 0));
        
        dlmwriteVec(inFilePath, outFilePath, numPlants, num_columns, gap, burn,
                    numSamples);
    }
  
}
