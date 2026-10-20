# OME basic plotting and analysis functions for sequencing data


#' Check available packages are installed and loaded
#' 
#' @description
#' Checks packages are installed and loaded
#' 
#' @param required_pkgs list of required packages.
#' 
#' @returns Loads and installs packages.
#' 
#' @export
#' 
#' @examples 
#' pkg_check(c("phyloseq"))
pkg_check <- function(required_pkgs = NULL) {
  
  # required packages for function
  all_pkgs <- c("tidyverse","here","phyloseq",required_pkgs)
  
  # check BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  library("BiocManager")
  
  # check packages required are downloaded
  for (pkg in all_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      sprintf("The '%s' package is required, installing now.", pkg)
      BiocManager::install(pkg, dependencies = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Check if microViz is installed
#' 
#' @description
#' Checks if microViz is installed, if not installs it, and loads library either way.
#' 
#' @returns Loads and installs MicroViz
#' 
#' @export
#' 
#' @examples 
#' microviz_check()
microviz_check <- function() {
  
  # check BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  library("BiocManager")
  
  # check microviz required packages are downloaded
  pkg_check(c("phyloseq", "microbiome", "ComplexHeatmap"))
  
  # checks if microviz is installed
  if (!requireNamespace("microViz", quietly = TRUE)) {
    install.packages(
    "microViz",
    repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos")))
  }
  library("microViz")
}

  
#' Load a taxa table from decontam output
#' 
#' @description
#' This function reads in a taxa table from decontam output.
#' 
#' @param tax_path Path to the Taxonomy table (CSV)
#' @param taxa_level default is null, if taxa_level = "nine" will split for pr2 levels.
#'  
#' @returns a taxa table object.
#' 
#' @export
#' 
#' @examples 
#' read_tax_table(here("MiFish_final_taxonomy_table.csv", taxa_level = NULL))
read_tax_table <- function(tax_path, taxa_level = NULL) {
  
  # required packages for function
  required_pkgs <- c("tidyverse","here")
  
  # required packages for function
  pkg_check()
  
  # load in taxa table
  # taxa table
  tax_table <- read_csv(here(tax_path)) 
  
  if(!is.null(taxa_level)) { 
    if(taxa_level == "nine") {
    tax_table <-tax_table %>%
    separate_wider_delim(sum.taxonomy, 
                         delim = ";", 
                         names = c("Domain", "Supergroup", "Division", "Subdivision", 
                                   "Class","Order", "Family", "Genus","Species"))
    } else { message("error: taxa_level either nine or null") }
  } else {
    tax_table <-tax_table %>%
      separate_wider_delim(sum.taxonomy, 
                           delim = ";", 
                           names = c("Domain", "Phylum","Class","Order", "Family", "Genus","Species"))
    
  }
  tax_table <- tax_table %>%
    column_to_rownames("seq_id") 
  
  return(tax_table)
}

#' Load asv table from decontam output
#' 
#' @description
#' This function reads in a asv table from decontam output.
#' 
#' @param asv_path Path to the asv table (csv)
#' 
#' @returns an asv table object.
#' 
#' @export
#' 
#' @examples 
#' read_asv_table(here("MiFish_final_asv_table.csv"))
read_asv_table <- function(asv_path) {
  
  # required packages for function
  pkg_check()
  
  # load in asv table
  otu_table <- read_csv(here(asv_path)) %>%
    column_to_rownames("seq_id")
  
  return(otu_table)
}

#' Load metadata from decontam output
#' 
#' @description
#' This function reads metadata from decontam output.
#' 
#' @param meta_path Path to the metadata
#' @param col_select columns to select for analysis
#' 
#' @returns metadata
#' 
#' @export
#' 
#' @examples 
#' meta <- read_metadata(here("MiFish_final_metadata.csv"))
read_metadata <- function(meta_path, col_select = NULL) {
  
  # required packages for function
  pkg_check()
  
  # load in asv table
  metadata <- read_csv(here(meta_path))
  
  # if col_select is not NULL, select columns
  if(!is.null(col_select)) {
    metadata <- metadata %>%
      dplyr::select(all_of(col_select))
  } else {
    message("No metadata columns selected")
  }
  
  return(metadata)
}

#' Load fasta file
#' 
#' @description
#' This function reads in a fasta file so sequences are matched with seq_id.
#' 
#' @param fasta_path Path to the fasta file.
#' 
#' @returns fasta seqs and seq_ids
#' 
#' @export
#' 
#' @examples 
#' read_fasta("MiFish_final_seqs.fasta")
read_fasta <- function(fasta_path) {
  
  # required packages for function
  pkg_check(c("Biostrings"))
  
  # load in asv table
  seqs <- readDNAStringSet(fasta_path)
  # format for use in R:
  seq_df <- data.frame(
    seq_id = names(seqs),
    sequence = as.character(seqs),
    seq_length = width(seqs),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  return(seq_df)
}

#' Load and build a phyloseq object
#' 
#' @description
#' This function reads in sequencing data (otu, tax, meta, fasta) and 
#' produces a phyloseq object.
#' 
#' @param asv_path Path to the OTU/ASV table (CSV)
#' @param tax_path Path to the Taxonomy table (CSV)
#' @param tax_level If this is = "nine" will correct for pr2 9 level, 
#' otherwise separates 7 level taxonomy
#' @param meta_path Path to the Metadata file (CSV)
#' @param fasta_path Path to fasta table (.fasta)
#' @param col_select Character vector of columns to keep from metadata.
#' @param output_dir Directory where the phyloseq.rds file will be saved
#' @param filename filename to save phyloseq rds object as.
#' 
#' @returns a phyloseq object, phyloseq object saved in an .rds file
#' 
#' @export
#' 
#' @examples 
#' create_phyloseq_object(asv_path = here("MiFish_final_asv_table.csv"),
#'                         tax_path = here("MiFish_final_taxonomy_table.csv"),
#'                         meta_path = here("MiFish_final_metadata.csv"),
#'                         fasta_path = here("MiFish_final_seqs.fasta"),
#'                         # if you want to choose particular columns for metadata 
#'                         # replace null with list of column id's e.g. c("samp_name", "station")
#'                         col_select = NULL,
#'                         output_dir = "",
#'                         filename = "MiFish_phyloseq_obj")
#'                         
create_phyloseq_object <- function(asv_path, tax_path, taxa_level = NULL,
                                   meta_path, fasta_path,
                                   col_select = NULL, output_dir, filename) {
  
  # required packages for function
  pkg_check("Biostrings")
  
  # load in sequencing data using above functions:
  # asv table
  asv_data <- read_asv_table(here(asv_path))
  # taxa table
  tax_data <- read_tax_table(here(tax_path), taxa_level = NULL)
  # metadata table
  metadata <- read_metadata(here(meta_path), col_select) %>% 
    column_to_rownames("samp_name")
  # fasta seqs
  fasta_seqs <- readDNAStringSet(fasta_path)
  
  # Format data for phyloseq:
  ASV <- otu_table(as.matrix(asv_data), taxa_are_rows = TRUE)
  TAX <- tax_table(as.matrix(tax_data))
  META <- sample_data(metadata) 
  SEQS <- refseq(fasta_seqs)
  
  # Create phyloseq object:
  ps_obj <- phyloseq(ASV, TAX, META, SEQS)
  
  # Create and save asv id for physeq and matching hashmd5 
  dummy_ids <- data.frame(hash_md5 = taxa_names(ps_obj),
                          ASV_dummy = paste0("ASV", "_",(1:ntaxa(ps_obj))))
  
  # Replace the taxa names across the entire phyloseq object
  taxa_names(ps_obj) <- dummy_ids$ASV_dummy
  # save key
  write_csv(dummy_ids, here(output_dir,paste0(filename,"_hashkey.csv")))
  
  # save phyloseq object for later use:
  saveRDS(ps_obj, here(output_dir,paste0(filename,".rds")))
  
  # output ps object:
  print(ps_obj)
}

#' Create figures directory
#' 
#' @param folder_name folder name for figure output if different than default. 
#' (default: figures)
#' 
#' @description
#' This function creates a directory to save figures.
#'
#' @examples 
#' create_dir()
create_dir <- function(folder_name = "figures") {
  
  # required packages for function
  pkg_check()
  
  # check if directory exists and if not, create it
  if (!dir.exists(here::here(folder_name))) {
    dir.create(here::here(folder_name), recursive = TRUE)
  }
}

