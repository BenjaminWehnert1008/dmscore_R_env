# input: wildtype-seq (string or fasta file), start&stopp pos., output_folder_path
# output: .csv with all possible NNK codons for each position


suppressMessages(library(Biostrings))

# Define the function
generate_possible_variants <- function(wt_seq_input, start_stop_pos, output_file) {
  # Parse the start and stop positions from the input format "start-stop"
  positions <- unlist(strsplit(start_stop_pos, "-"))
  start_pos <- as.numeric(positions[1])
  stop_pos <- as.numeric(positions[2])

  # Check if the input is a file or a string
  if (file.exists(wt_seq_input)) {
    # If it's a file, read the sequence from the fasta file
    seq_data <- readDNAStringSet(filepath = wt_seq_input)
    wt_seq <- seq_data[[1]]  # Extract the sequence
  } else {
    # Otherwise, treat the input as a sequence string
    wt_seq <- DNAString(wt_seq_input)
  }


  # Extract the sequence between start and stop codons
  coding_seq <- subseq(wt_seq, start = start_pos, end = stop_pos)

  class(coding_seq)

  # List of predefined NNK codons
  nnk_codons <- c('AAG', 'AAT', 'ATG', 'ATT', 'AGG', 'AGT', 'ACG', 'ACT',
                  'TAG', 'TAT', 'TTG', 'TTT', 'TGG', 'TGT', 'TCG', 'TCT',
                  'GAG', 'GAT', 'GTG', 'GTT', 'GGG', 'GGT', 'GCG', 'GCT',
                  'CAG', 'CAT', 'CTG', 'CTT', 'CGG', 'CGT', 'CCG', 'CCT')

  coding_seq <- as.character(coding_seq)

  # Function to split a DNA sequence into codons (triplets)
  split_into_codons <- function(seq) {
    return(strsplit(seq, "(?<=.{3})", perl = TRUE)[[1]])
  }

  # Split wild-type sequence into codons
  wt_codons <- split_into_codons(coding_seq)

  # Initialize DataFrame to store mutated variants
  result <- data.frame(Codon_Number = integer(), wt_codon = character(), Variant = character(), stringsAsFactors = FALSE)

  # Iterate over each codon in the wild-type sequence
  for (i in seq_along(wt_codons)) {
    wt_codon <- wt_codons[i]

    # Filter NNK codons that are different from the wild-type codon
    possible_variants <- nnk_codons[nnk_codons != wt_codon]

    # Add all variants to the result list, including the wild-type codon
    for (variant in possible_variants) {
      result <- rbind(result, data.frame(Codon_Number = i, wt_codon = wt_codon, Variant = variant, stringsAsFactors = FALSE))
    }
  }

  # Save the variants into a CSV file
  write.csv(result, output_file, row.names = FALSE)

}

# example
# generate_possible_variants("/Users/benjaminwehnert/CRG/DMS_QC/testing_data/MORtn5_reference.fa", "23-1225", "/Users/benjaminwehnert/CRG/DMS_QC/testing_data/testing_outputs/possible_NNK_mutations.csv")













