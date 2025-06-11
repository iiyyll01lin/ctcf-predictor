# Load aligned sequences (use smaller subset to avoid memory issues)
cat("Loading aligned sequences...\n")

# Check if we have a small subset available, otherwise create one
if (!file.exists('data/aligned_subset_small.fasta')) {
  cat("Creating smaller subset for testing...\n")
  system("head -n 2000 data/aligned_sequences.fasta > data/aligned_subset_small.fasta")
}

# Use the smaller subset to avoid hanging
aligned_seqs <- readDNAStringSet('data/aligned_subset_small.fasta')
cat("Loaded", length(aligned_seqs), "aligned sequences from subset\n")