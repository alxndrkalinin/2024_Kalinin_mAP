suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(readr))

# Function to initialize the Seurat object from 10X data
initialize_seurat <- function(gse_id, base_dir) {
  data_dir <- file.path(base_dir, gse_id)
  cat("Loading 10X data from:", data_dir, "\n")
  perturbseq_data <- Seurat::Read10X(data.dir = data_dir)
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = perturbseq_data,
    project = "crispri",
    min.cells = 3,
    min.features = 200
  )
  cat("Created Seurat object with", ncol(seurat_obj), "cells\n")
  
  # Calculate mitochondrial percentage
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  return(seurat_obj)
}

# Function to filter cells based on RNA counts and mitochondrial percentage
filter_cells <- function(seurat_obj) {
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > 2000 & nFeature_RNA < 7000 & percent.mt < 15
  )
  cat("After filtering, object has", ncol(seurat_obj), "cells\n")
  return(seurat_obj)
}

# Function to normalize and scale the data
normalize_and_scale <- function(seurat_obj) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  all.genes <- rownames(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = all.genes)
  return(seurat_obj)
}

# Function to identify variable features (without plotting)
identify_variable_features <- function(seurat_obj, nfeatures = 1000, top_n = 10) {
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures)
  top10 <- head(VariableFeatures(seurat_obj), top_n)
  cat("Top", top_n, "variable features:\n")
  print(top10)
  return(seurat_obj)
}

# Function to extract scaled data for selected variable features and save to file
extract_and_save_scaled_data <- function(seurat_obj, processed_output_file) {
  scaled_data <- Seurat::GetAssayData(object = seurat_obj, slot = "scale.data")
  selected_features <- VariableFeatures(seurat_obj)
  
  scaled_df <- as.data.frame(scaled_data) %>%
    tibble::rownames_to_column(var = "gene") %>%
    dplyr::filter(gene %in% selected_features)
  
  cat("Scaled data dimensions:", dim(scaled_df), "\n")
  
  readr::write_tsv(scaled_df, processed_output_file)
  cat("Saved processed data to", processed_output_file, "\n")
}

# Main function to run the data processing pipeline (without plotting)
main <- function() {
  gse_id <- "GSE132080"
  base_dir <- file.path("inputs")
  output_dir <- "outputs"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  processed_output_file <- file.path(output_dir, paste0(gse_id, "_processed_matrix.tsv.gz"))
  cat("Output file:", processed_output_file, "\n")
  
  # Initialize Seurat object from 10X data
  seurat_obj <- initialize_seurat(gse_id, base_dir)
  print(seurat_obj)
  
  # Filter cells
  seurat_obj <- filter_cells(seurat_obj)
  
  # Normalize and scale the data
  seurat_obj <- normalize_and_scale(seurat_obj)
  
  # Identify variable features (without plotting)
  seurat_obj <- identify_variable_features(seurat_obj, nfeatures = 1000, top_n = 10)
  
  # Extract and save the scaled data for the selected features
  extract_and_save_scaled_data(seurat_obj, processed_output_file)
}

main()