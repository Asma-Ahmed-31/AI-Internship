# =====================================================================
#               AI and Biotechnology / Bioinformatics
# =====================================================================

# ---------------------------------------------------------------------
#              AI and Omics Data Analysis (Microarray)
# ---------------------------------------------------------------------

# 0. Install and Load Required Packages
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install(c("GEOquery","affy","arrayQualityMetrics"), ask = FALSE)
install.packages("dplyr")

library(GEOquery)
library(affy)
library(arrayQualityMetrics)
library(dplyr)

# -------------------------------------
# Download Series Matrix Files
# -------------------------------------
gse_list <- getGEO("GSE79973", GSEMatrix = TRUE)
# getGEO returns a list; select the first ExpressionSet
gse_data <- if (is.list(gse_list)) gse_list[[1]] else gse_list

# Extract expression / feature / phenotype
expression_data <- exprs(gse_data)
feature_data    <- fData(gse_data)
phenotype_data  <- pData(gse_data)

# Check missing values in sample annotation
sum(is.na(phenotype_data$source_name_ch1))

# --------------------------------------
# Download Raw Data (CEL files)
# --------------------------------------
# Attempt to download supplements but don't stop script if it fails
tryCatch({
  getGEOSuppFiles("GSE79973", baseDir = "Raw_Data", makeDirectory = TRUE)
}, error = function(e) {
  message("Warning: getGEOSuppFiles() failed or no supplements available. Continuing with series matrix only.")
})

# Typical expected archive paths
tar_path <- file.path("Raw_Data", "GSE79973", "GSE79973_RAW.tar")
zip_path <- file.path("Raw_Data", "GSE79973", "E-GEOD-79973.zip")

# If tar exists, untar; if zip exists, unzip. Only proceed if extraction succeeded.
if (file.exists(tar_path)) {
  untar(tar_path, exdir = "Raw_Data/CEL_Files")
} else if (file.exists(zip_path)) {
  unzip(zip_path, exdir = "Raw_Data/E_GEOD79973")
} else {
  # also check if a tar or zip is directly in Raw_Data folder
  alt_tar <- list.files("Raw_Data", pattern = "\\.tar(\\.gz)?$", full.names = TRUE)
  alt_zip <- list.files("Raw_Data", pattern = "\\.zip$", full.names = TRUE)
  if (length(alt_tar) > 0) untar(alt_tar[1], exdir = "Raw_Data/CEL_Files")
  if (length(alt_zip) > 0) unzip(alt_zip[1], exdir = "Raw_Data/E_GEOD79973")
}

# Read CEL files only if present
cel_files <- list.files("Raw_Data/CEL_Files", pattern = "\\.[cC][eE][lL]$", recursive = TRUE, full.names = TRUE)
if (length(cel_files) > 0) {
  raw_data <- ReadAffy(celfile.path = "Raw_Data/CEL_Files")
  raw_data   # Displays basic information about the dataset
  
  # QC on raw data (wrap in tryCatch to avoid stopping on QC issues)
  tryCatch({
    arrayQualityMetrics(expressionset = raw_data,
                        outdir = "Results/QC_Raw_Data",
                        force = TRUE,
                        do.logtransform = TRUE)
  }, error = function(e) message("arrayQualityMetrics (raw) failed: ", e$message))
  
  # RMA normalization on raw data
  normalized_data <- rma(raw_data)
  
  # QC after normalization
  tryCatch({
    arrayQualityMetrics(expressionset = normalized_data,
                        outdir = "Results/QC_Normalized_Data",
                        force = TRUE)
  }, error = function(e) message("arrayQualityMetrics (normalized) failed: ", e$message))
  
  processed_data <- as.data.frame(exprs(normalized_data))
} else {
  # Fall back to series-matrix expression data if no raw CELs available
  message("No CEL files found; using series-matrix expression data for downstream steps.")
  processed_data <- as.data.frame(expression_data)
}

dim(processed_data)   # Dimensions: number of probes × number of samples

# ---------------------------------------------------------------------------
# Filter Low-Variance / Low-Intensity Transcripts (“soft” intensity based)
# ---------------------------------------------------------------------------

# Use apply() to compute median per row (no extra packages required)
row_median <- apply(as.matrix(processed_data), 1, median, na.rm = TRUE)

# Visualize distribution of probe median intensities (optional)
hist(row_median, breaks = 100, freq = FALSE, main = "Median Intensity Distribution")
threshold <- 3.5
abline(v = threshold, col = "black", lwd = 2)

indx <- row_median > threshold
filtered_data <- processed_data[indx, ]

# Rename columns to phenotype sample names only if counts match
if (ncol(filtered_data) == nrow(phenotype_data)) {
  colnames(filtered_data) <- rownames(phenotype_data)
}

processed_data <- filtered_data

# -----------------------------------
# Phenotype Data Preparation
# -----------------------------------
class(phenotype_data$source_name_ch1)

groups <- factor(phenotype_data$source_name_ch1,
                 levels = c("gastric mucosa", "gastric adenocarcinoma"),
                 label = c("normal", "cancer"))

class(groups)
levels(groups)

# End of script
dim(filtered_data)


ls()

library(Biobase)
# Convert to ExpressionSet
eset_filtered <- ExpressionSet(
  assayData = as.matrix(filtered_data),
  phenoData = AnnotatedDataFrame(phenotype_data)
)

library(arrayQualityMetrics)

arrayQualityMetrics(
  expressionset = eset_filtered,
  outdir = "Quality_Control_Report",
  force = TRUE,
  do.logtransform = FALSE
)

arrayQualityMetrics(expressionset = eset_filtered,
                    outdir = "Quality_Control_Report",
                    force = TRUE,
                    do.logtransform = FALSE)




