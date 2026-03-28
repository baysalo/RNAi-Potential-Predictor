# ==============================================================================
# RNAi Potential Predictor: Bacterial-Fungal Cross-Kingdom Silencing
# Author: Ömür BAYSAL Ph.D. 
# License: MIT
# ==============================================================================

library(Biostrings)

# --- 1. SETTINGS & PATHS ---
bact_file <- "data/bacterial_genome.fasta"
fung_file <- "data/fungal_transcriptome.fasta"
host_file <- "data/host_genome.fasta"

# --- 2. LOADING DATA ---
message("Loading genomic sequences...")
bact_seqs  <- readDNAStringSet(bact_file)
fung_trans <- readDNAStringSet(fung_file)
host_gen   <- readDNAStringSet(host_file)

# --- 3. DICING FUNCTION ---
generate_siRNAs <- function(seq, size = 21) {
  if (nchar(seq) < size) return(NULL)
  views <- Views(seq, start = 1:(nchar(seq) - size + 1), width = size)
  return(as(views, "DNAStringSet"))
}

# --- 4. SCREENING LOOP (Pathogen Match) ---
all_results <- list()
message("Scanning bacterial genome for fungal matches...")

for(i in seq_along(bact_seqs)) {
  current_siRNAs <- generate_siRNAs(bact_seqs[[i]])
  if(is.null(current_siRNAs)) next
  
  # Map to fungus (0 mismatches for high potency)
  hits <- vcountPDict(current_siRNAs, fung_trans, max.mismatch = 0)
  valid_idx <- which(rowSums(hits) > 0)
  
  if(length(valid_idx) > 0) {
    all_results[[i]] <- current_siRNAs[valid_idx]
  }
}

effective_siRNAs <- unique(do.call(c, all_results))

# --- 5. HOST SAFETY FILTERING ---
message("Filtering for host off-targets (Safety Check)...")
host_hits <- vcountPDict(effective_siRNAs, host_gen, max.mismatch = 2)
safe_indices <- which(rowSums(host_hits) == 0)

if(length(safe_indices) > 0) {
  final_safe_siRNAs <- effective_siRNAs[safe_indices]
} else {
  stop("No safe candidates found.")
}

# --- 6. TARGET MAPPING & EXPORT ---
message("Finalizing target identification...")
hits_matrix <- vcountPDict(final_safe_siRNAs, fung_trans, max.mismatch = 0)
hit_pairs <- which(hits_matrix > 0, arr.ind = TRUE)

final_mapping <- data.frame(
  siRNA_ID = paste0("siRNA_", hit_pairs[,1]),
  siRNA_Sequence = as.character(final_safe_siRNAs[hit_pairs[,1]]),
  Target_Fungal_Transcript = names(fung_trans)[hit_pairs[,2]],
  GC_Percent = as.numeric(letterFrequency(final_safe_siRNAs[hit_pairs[,1]], "GC", as.prob = TRUE)) * 100
)

# Export Results
write.csv(final_mapping, "results/Final_RNAi_Mapping_Results.csv", row.names = FALSE)
writeXStringSet(final_safe_siRNAs, "results/Final_Safe_siRNAs.fasta")

message("Pipeline Complete. Check the /results folder.")
