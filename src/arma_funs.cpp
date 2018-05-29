
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <math.h> 
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

// 
// [[Rcpp::export]]
Rcpp::IntegerVector countNonzeroColsArma (arma::mat mat) {
  
  int n = mat.n_rows, m = mat.n_cols;
  
  Rcpp::IntegerVector nonzeroCount(m, 0);
  
  // initialize progress bar
  Progress p(m, true);
  
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      if (mat(i, j) > 0) {
        nonzeroCount[j]++;
      }
    } // end loop through columns
    // update progress
    p.increment(); 
  } // end loop through rows
  
  return(nonzeroCount);
  
}

// [[Rcpp::export]]
IntegerVector countNonzeroCols (NumericMatrix mat) {
  
  int n = mat.nrow(), m = mat.ncol();
  
  IntegerVector nonzeroCount(m, 0);
  
  // initialize progress bar
  Progress p(m, true);
  
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      if (mat(i, j) > 0) {
        nonzeroCount[j]++;
      }
    } // end loop through columns
    // update progress
    p.increment(); 
  } // end loop through rows
  
  return(nonzeroCount);
  
}

/* 
* Column medians 
*/

// [[Rcpp::export]]
NumericVector colMediansC(NumericMatrix x) {
  int ncol = x.ncol(); 
  NumericVector out(ncol);
  
  NumericVector col_j;
  
  for (int j = 0; j < ncol; j++) {
    col_j = x(_, j);
    out[j] = median(col_j);
  }
  
  return out;
  
}

/*  //---------------------------------------------------------------------
* TIC normalization
*/ //---------------------------------------------------------------------

// [[Rcpp::export]]
NumericVector rowMeansC(NumericMatrix x) {
  
  return rowMeans(x);
  
}

// [[Rcpp::export]]
List normalizeTIC(NumericMatrix &dataMatrix) {
  
  // Copy matrix without propogating changes to original for the output
  int nrow = dataMatrix.nrow();
  int ncol = dataMatrix.ncol();
  NumericMatrix dataNormPixel = dataMatrix(Range(0, nrow - 1), 
                                           Range(0, ncol - 1));
  
  // Normalization factor is the row means
  NumericVector normFactorPixel = rowMeans(dataMatrix);
  
  // Divide each row of the datamatrix by its row mean
  for (int i = 0; i < nrow; i++) {
    dataNormPixel(i, _) = dataNormPixel(i, _) / normFactorPixel[i];
  }
  
  // Give output
  return List::create(
    _["dataNormPixel"] = dataNormPixel,
    _["normFactorPixel"] = normFactorPixel
  );
  
}



/*  //---------------------------------------------------------------------
* Median normalization
*/ //---------------------------------------------------------------------

// [[Rcpp::export]]
double nonzeroMedianC(NumericVector x) {
  NumericVector keep = x[x != 0];
  
  return median(keep);
  
}

// [[Rcpp::export]]
NumericVector rowNonzeroMediansC(NumericMatrix &x) {
  int nrow = x.nrow(); 
  NumericVector out(nrow);
  
  NumericVector row_i;
  
  for (int i = 0; i < nrow; i++) {
    row_i = x(i, _);
    out[i] = nonzeroMedianC(row_i);
  }
  
  return out;
  
}

// [[Rcpp::export]]
List normalizeNonzeroMedian(NumericMatrix &dataMatrix) {
  
  // Copy matrix without propogating changes to original
  int nrow = dataMatrix.nrow();
  int ncol = dataMatrix.ncol();
  NumericMatrix dataNormPixel = dataMatrix(Range(0, nrow - 1), 
                                           Range(0, ncol - 1));
  
  // Calculate normalization factors
  NumericVector normFactorPixel = rowNonzeroMediansC(dataMatrix);
  
  // Loop through rows and divide by nonzero median
  for (int i = 0; i < nrow; i++) {
    dataNormPixel(i, _) = dataNormPixel(i, _) / normFactorPixel[i];
  }
  
  return List::create(
    _["dataNormPixel"] = dataNormPixel,
    _["normFactorPixel"] = normFactorPixel
  );
  
  
  
  //
  
}


/*  //---------------------------------------------------------------------
* Median fold change
*/ //---------------------------------------------------------------------


// [[Rcpp::export]]
NumericVector rowMFCsC (NumericMatrix dataMatrix, 
                        NumericVector refSpec) {
  // Calculate median fold change between every row of data matrix 
  // and a reference spectrum. Only median fold change between elements of 
  // refSpec which are nonzero
  
  // Number of rows in data matrix
  int nrow = dataMatrix.nrow();
  
  // Which indices of refSpec are nonzero, and retain only those values
  LogicalVector keepBool = (refSpec > 0);
  IntegerVector keepIdx = seq(0, refSpec.length() - 1);
  keepIdx = keepIdx[keepBool];
  
  // Initialize output, a vector of median fold changes
  NumericVector MFCvec(nrow);
  
  // Initialize a vector for each row 
  NumericVector row_i;
  
  // Calculate
  for (int i = 0; i < nrow; i++) {
    
    // i-th row
    row_i = dataMatrix(i, _);
    
    // Median fold change
    MFCvec[i] = median(row_i[keepIdx] / refSpec[keepIdx]);
  }
  
  // Output
  return MFCvec;
  
}


// [[Rcpp::export]]
NumericVector rowMFCsC_nonzero (NumericMatrix dataMatrix, 
                                NumericVector refSpec) {
  // Calculate median of *nonzero* fold change between every row of data matrix 
  // and a reference spectrum. Only median fold change between elements of 
  // refSpec which are nonzero
  
  // Number of rows in data matrix (spectra)
  int nrow = dataMatrix.nrow();
  
  // Which indices of reference spectrum
  LogicalVector keepBoolRefSpec = (refSpec > 0);
  
  // Initialize output
  NumericVector MFCvec(nrow);
  
  // Initialize a vector for each row 
  NumericVector row_i;
  
  // Initialize vectors for which indices to consider for each row
  LogicalVector keepBoolRow;
  IntegerVector KeepIdxRow;
  
  for (int i = 0; i < nrow; i++) {
    
    // Only keep elements of the row which are nonzero
    keepBoolRow = (dataMatrix(i, _) > 0);
    KeepIdxRow = seq(0, refSpec.length() - 1);
    KeepIdxRow = KeepIdxRow[keepBoolRow & keepBoolRefSpec];
    
    // i-th row
    row_i = dataMatrix(i, _);
    
    // median of the nonzero foldchanges
    MFCvec[i] = median(row_i[KeepIdxRow] / refSpec[KeepIdxRow]);
    
  }
  
  return MFCvec;
  
}



// // [[Rcpp::export]]
// NumericMatrix normalizeMFC(NumericMatrix dataMatrix, NumericVector ref) {
//   // use median or non-zero median fold change?
//  
//   int nrow = x.nrow(), ncol = x.ncol();
//   NumericMatrix out = x(Range(0, nrow - 1), Range(0, ncol - 1));
//  
//  
//   // check which mz's in reference to use for 
//   LogicalVector keep_logical = (ref > 0);
//   IntegerVector keep_index = seq(0, ref.length() - 1);
//   keep_index = keep_index[keep_logical];
// 
//   
//   NumericVector factors = rowMFCsC(x, ref);
//   
//   for (int i = 0; i < nrow; i++) {
//     out(i, _) = x(i, _) / factors[i];
//   }
//   
//   return out;
//   
// }





// [[Rcpp::export]]
List normalizeMFC_nonzero(NumericMatrix dataMatrix, NumericVector refSpec) {
  // Median fold-change normalization
  
  // Copy matrix without propogating changes to original
  int nrow = dataMatrix.nrow(), ncol = dataMatrix.ncol();
  NumericMatrix dataNormPixel = dataMatrix(Range(0, nrow - 1), 
                                           Range(0, ncol - 1));
  
  
  // Normalization
  NumericVector normFactorPixel = rowMFCsC_nonzero(dataMatrix, refSpec);
  
  for (int i = 0; i < nrow; i++) {
    dataNormPixel(i, _) = dataNormPixel(i, _) / normFactorPixel[i];
  }
  
  return List::create(
    _["dataNormPixel"] = dataNormPixel,
    _["normFactorPixel"] = normFactorPixel
  );
}

/*  //---------------------------------------------------------------------
*  Inter-sample normalization 
*/ // --------------------------------------------------------------------

// [[Rcpp::export]]
List intersampleNormTIC(NumericMatrix dataNormPixel, 
                        IntegerVector sampleFirstIndices) {
  
  // sampleFirstIndices is a vector of the first index of the data matrix
  // corresponding to each sample (patient)
  
  // Allocate output (normalized data)
  int nrow = dataNormPixel.nrow(), ncol = dataNormPixel.ncol();
  NumericMatrix dataNormSample(nrow, ncol);
  
  // Number of samples
  int m = sampleFirstIndices.length();
  
  // Allocate vector for per-sample normalization factors
  NumericVector normFactorSample(m);
  
  // Allocate space for data_matrix and sample profile spectrum for each sample
  NumericMatrix dataMatrixK;
  NumericVector sampleProfile_k;
  
  // Allocate space for first and last indices for each patient
  int first, last;
  
  // loop through all the patients
  for (int j = 0; j < m; j++) {
    
    // Grab row indices for this patient
    first = sampleFirstIndices[j] - 1;
    // Check if we're at the last patient
    if (j < m - 1) { 
      last  = sampleFirstIndices[j + 1] - 1;
    } else {
      last  = dataNormPixel.nrow() - 1;
    }
    
    // Data submatrix for this patient
    dataMatrixK = dataNormPixel(Range(first, last), _ ); 
    
    // Calculate representative spectrum for this sample
    sampleProfile_k = colMeans(dataMatrixK);
    
    // Caclulate normalization factor
    normFactorSample[j] = mean(sampleProfile_k);
    
    // Fill out the output matrix by looping through rows of patient data
    for (int row = 0; row <= (last - first); row++) {
      dataNormSample(first + row, _) = 
        dataMatrixK(row, _) / normFactorSample[j];
    }
    
  }
  
  return List::create(
    _["dataNormSample"] = dataNormSample,
    _["normFactorSample"] = normFactorSample
  );
  
  
}

// [[Rcpp::export]]
List intersampleNormMedian(NumericMatrix dataNormPixel, 
                           IntegerVector sampleFirstIndices) {
  
  // sampleFirstIndices is a vector of the first index of the data matrix
  // corresponding to each sample (patient)
  
  // Allocate output (normalized data)
  int nrow = dataNormPixel.nrow(), ncol = dataNormPixel.ncol();
  NumericMatrix dataNormSample(nrow, ncol);
  
  // Number of samples
  int m = sampleFirstIndices.length();
  
  // Allocate vector for per-sample normalization factors
  NumericVector normFactorSample(m);
  
  // Allocate space for data_matrix and reference spectrum for each sample
  NumericMatrix dataMatrixK;
  NumericVector sampleProfile_k;
  
  // IntegerVector sampleRowIndices;
  int first, last;
  
  // loop through all the patients
  for (int j = 0; j < m; j++) {
    
    // Grab row indices for this patient
    first = sampleFirstIndices[j] - 1;
    // Check if we're at the last patient
    if (j < m - 1) { 
      last  = sampleFirstIndices[j + 1] - 1;
    } else {
      last  = dataNormPixel.nrow() - 1;
    }
    
    // Data submatrix for this patient
    dataMatrixK = dataNormPixel(Range(first, last), _ ); 
    
    // Calculate profile spectrum for this sample
    sampleProfile_k = colMeans(dataMatrixK);
    
    // Caclulate normalization factor (median)
    normFactorSample[j] = nonzeroMedianC(sampleProfile_k);
    
    // Fill out the output matrix by looping through rows of patient data
    for (int row = 0; row <= (last - first); row++) {
      dataNormSample(first + row, _) = dataMatrixK(row, _) / normFactorSample[j];
    }
    
  }
  
  return List::create(
    _["dataNormSample"] = dataNormSample,
    _["normFactorSample"] = normFactorSample
  );
  
  
  
}

// [[Rcpp::export]]
List intersampleNormMFC(NumericMatrix x, IntegerVector sampleFirstIndices) {
  
  // sampleFirstIndices is a vector of the first index of the data matrix 
  //      (in R index mode, so subtract 1!)
  // corresponding to each sample (patient)
  
  // Allocate output (normalized data)
  int nrow = x.nrow(), ncol = x.ncol();
  NumericMatrix dataNormSample(nrow, ncol);
  
  // Number of samples
  int m = sampleFirstIndices.length();
  
  // Allocate vector for per-sample normalization factors
  NumericVector normFactorSample(m);
  
  // Allocate space for data_matrix and reference spectrum for each sample
  NumericMatrix dataMatrixK;
  NumericMatrix sampleProfiles(m, ncol);
  
  // Vectors for first and last indices for each patient
  IntegerVector firstVec(m), lastVec(m);
  
  // loop through all the patients and create sample profile
  for (int j = 0; j < m; j++) {
    
    // Grab row indices for this patient
    firstVec[j] = sampleFirstIndices[j] - 1;
    // Check if we're at the last patient
    if (j < m - 1) { 
      lastVec[j]  = sampleFirstIndices[j + 1] - 1;
    } else {
      lastVec[j]  = x.nrow() - 1;
    }
    
    // Data submatrix for this patient
    dataMatrixK = x(Range(firstVec[j], lastVec[j]), _ ); 
    
    // Calculate representative spectrum for this sample
    sampleProfiles(j, _) = colMeans(dataMatrixK);
    
  }
  
  // Reference spectrum for all samples; mean 
  NumericVector refSpecAll = colMediansC(sampleProfiles);
  
  // Normalization factor (nonzero median fold change)
  normFactorSample = rowMFCsC_nonzero(sampleProfiles, refSpecAll);
  
  // first and last index for a specific patient 
  int firstEll, lastEll;
  
  // Loop through patients again
  for (int ell = 0; ell < m; ell++) {
    
    firstEll = firstVec[ell];
    lastEll = lastVec[ell];
    
    // Fill out the output matrix by looping through rows of patient data
    for (int row = firstEll; row <= lastEll; row++) {
      dataNormSample(row, _) = x(row, _) / normFactorSample[ell];
    }
    
    
  }
  
  // Output
  return List::create(
    _["dataNormSample"] = dataNormSample,
    _["normFactorSample"] = normFactorSample
  );
  
  
}



