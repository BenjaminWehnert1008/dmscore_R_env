# Input: pre_processed_raw_gatk_path, mutation library used (NNK/NNN/custom_codon_vector), output_path
# Output: gatk table filtered for only single-codon mutations that are part of the library

# Load necessary library
library("dplyr")
library("stringr")

filter_gatk_by_codon_library <- function(gatk_file_path, codon_library, output_file_path) {
  # Load the GATK table from the provided file path
  gatk_table <- read.csv(gatk_file_path)

  # Predefined codon libraries
  nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT',
                  'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT',
                  'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT',
                  'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')

  nnn_codons <- c('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT',
                  'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT',
                  'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT',
                  'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT',
                  'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT',
                  'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT',
                  'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT',
                  'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT')

  # Check if the codon_library is a predefined string or a custom vector
  if (is.character(codon_library) && length(codon_library) == 1) {
    if (codon_library == "NNK") {
      codon_set <- nnk_codons
    } else if (codon_library == "NNN") {
      codon_set <- nnn_codons
    } else {
      stop("Invalid predefined codon library specified.")
    }
  } else if (is.vector(codon_library)) {
    codon_set <- codon_library
  } else {
    stop("Invalid codon library format. Must be a predefined string or a custom vector.")
  }

  # Filter for single-codon mutations that are part of the codon library and make sure to handle some mistaken formatting from gatk (there are cases where more than 3 Bases are mutated, but column varying_codons == 1)
  filtered_gatk <- gatk_table %>%
    filter(varying_codons == 1 & sub(".*>", "", codon_mut) %in% codon_set) %>%
    rowwise() %>%
    filter({
      # Split base_mut into individual mutations
      mutations <- unlist(strsplit(base_mut, ",\\s*"))  # Splits by comma and removes extra spaces
      # Extract numeric positions from each mutation string
      positions <- as.numeric(str_extract(mutations, "^[0-9]+"))

      # Calculate the distance between the first and last position
      distance <- max(positions, na.rm = TRUE) - min(positions, na.rm = TRUE)

      # Keep rows where the distance is <= 2
      distance <= 2
    }) %>%
    ungroup()


  # Write the filtered GATK table to the output file path
  write.csv(filtered_gatk, file = output_file_path, row.names = FALSE)
}

# Example usage with predefined library (NNK):
# filtered_gatk <- filter_gatk_by_codon_library(raw_gatk, codon_library = "NNK")

# Example usage with custom library:
# custom_library <- c('AAG', 'TGT', 'GCT')
# filtered_gatk <- filter_gatk_by_codon_library(raw_gatk, custom_library = custom_library)

#filter_gatk_by_codon_library("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/raw_gatk.csv", codon_library = "NNK", output_file_path = "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/gatk_filtered_by_codon_library.csv")


