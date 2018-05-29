


require(dplyr)


normalizeData <- function(dataMatrix, 
                          pixelInfo,
                          mzVec = NA, 
                          perPixelMethod = c("TIC", "median", "MFC"),
                          perSampleMethod = c("TIC", "median", "MFC"),
                          refType = c("median", "mean"),
                          useLog = F,
                          startedLog = F,
                          percNonzeroCutoff = 0.10,
                          restrictLower = NA,
                          restrictUpper = NA) {
  
  # INPUTS: 
  #   dataMatrix - the feature matrix
  #   pixelInfo - data frame of pixel metadata
  #   mzVec - mz features, required only if restricting mass range
  #   perPixelMethod - pixel-specific normalization mode; 
  # 
  # OUTPUTS
  #   data_norm - final result
  
  # Standardize inputs by capitalizing them
  perPixelMethod <- toupper(perPixelMethod)
  perSampleMethod <- toupper(perPixelMethod)
  
  # Convert data to dense matrix
  dataMatrix <- as.matrix(dataMatrix)
  
  # Remove features
  countNonzeroCols(dataMatrix)
  
  # Check if restricted range requested, and apply restriction
  if (any(!is.na(c(restrictLower, restrictUpper)))) {
    
    if (any(is.na(mzVec))) {
      stop("Must provide mz vector if restricting mz range.")
    }
    
    if (!is.na(restrictLower) & is.na(restrictUpper)) {
      # Case 1: Lower limit only
      keep_features <- (mzVec > restrictLower)
      
    } else if (is.na(restrictLower) & !is.na(restrictUpper)) {
      # Case 2: Upper limit only
      keep_features <- (mzVec < restrictLower)
      
    } else {
      # Case 3: Both limits specified
      keep_features <- (mzVec > restrictLower & 
                          mzVec < restrictUpper)
      
    }
    
    dataMatrix <- dataMatrix[, keep_features]
    keep_mz <- mzVec[keep_features]
    
  } else {
    
    if (any(is.na(keep_mz))) {
      keep_mz <- "(not provided in function call)"
    } else {
      keep_mz <- mzVec
    }
  }
  
  
  # Per-pixel normalization
  perPixelNorm <- normalizeDataPixel(dataMatrix, perPixelMethod, refType)
  
  
  # Per-sample normalization, if requested
  if (length(perSampleMethod) == 1) {
    # Check if value provided is valid
    if (!(perSampleMethod %in% c("TIC", "MEDIAN", "MFC"))) {
      stop("Must provide valid per-sample normalization")
    }
    
    # Calculate first indices in data for each sample
    sampleFirstIndices <- pixelInfo %>% 
      group_by(k) %>% 
      summarize(n = n()) %>% 
      ungroup() %>% 
      mutate(cumsum = cumsum(n),
             lag_cumsum = lag(cumsum),
             first = ifelse(k == 1, 1L, lag_cumsum + 1L)) %>% 
      pull(first)
    
    if (perSampleMethod == "TIC") {
      perSampleNorm <- intersampleNormTIC(perPixelNorm$dataNormPixel, 
                                          sampleFirstIndices)
    } else if (perSampleMethod == "MEDIAN") {
      perSampleNorm <- intersampleNormMedian(perPixelNorm$dataNormPixel, 
                                             sampleFirstIndices)
    } else if (perSampleMethod == "MFC") { # Note mean of profile is used
      perSampleNorm <- intersampleNormMFC(perPixelNorm$dataNormPixel, 
                                          sampleFirstIndices)
    } else {
      stop("Must pick either TIC, median, of MFC per-sample normalization")
    }
    
    data_norm <- perSampleNorm$dataNormSample
    normFactorSample <- perSampleNorm$normFactorSample
    
  } else { # Case per
    data_norm <- perPixelNorm$dataNormPixel
    normFactorSample <- "("
  }
  
  # Take log, if requested
  if (useLog) {
    if (startedLog) {
      delta <- quantile(data_norm[data_norm != 0], 0.05) %>% as.mumeric()
    } else {
      delta <- 1
    }
    
    data_norm <- log(data_norm + 1)
    
  }
  
  # 
  return(list(
    "data_norm" = data_norm,
    "mzVec" = keep_mz,
    "normFactorPixel" = perPixelNorm$normFactorPixel,
    "normFactorSample" = normFactorSample
  ))
  
}

# Per-pixel normalization
normalizeDataPixel <- function(dataMatrix,
                               perPixelMethod = c("TIC", "median", "MFC"), 
                               refType = c("median", "mean")) {
  
  # Standardize inputs by capitalizing them
  perPixelMethod <- toupper(perPixelMethod)
  
  # Perform normalization
  if (perPixelMethod == "TIC") {
    
    dataNormPixel <- normalizeTIC(dataMatrix)
    
  } else if (perPixelMethod == "MEDIAN") {
    
    dataNormPixel <- normalizeNonzeroMedian(dataMatrix)
    
  } else if (perPixelMethod == "MFC") {
    
    if (length(refType) != 1) {
      stop("Must specify ONE OF mean or median reference for MFC norm!")
    }
    
    # Capitalize refType
    refType <- toupper(refType)
    
    # Check type of reference
    if (refType == "MEAN") {
      refSpec <- colMeansC(dataMatrix)
    } else if (refType == "MEDIAN") {
      refSpec <- colMediansC(dataMatrix)
    } else {
      stop("Must specify mean or median reference for MFC norm!")
    }
    
    # data_norm <- normalizeMFC(dataMatrix, ref) # Too many zeros...
    dataNormPixel <- normalizeMFC_nonzero(dataMatrix, refSpec)
    
  } else {
    
    stop("Error: must specify either TIC, median, or median fold change norm!")
    
  }
  
  
  # Output
  return(list(
    "dataNormPixel" = dataNormPixel$dataNormPixel,
    "normFactorPixel" = dataNormPixel$normFactorPixel
  ))
  
}


