#' Extract nucleotide contexts for each site in SNV hotspots.
#' 
#' @param seq Dataframe containing DNA sequences of each site in the SNV hotspot.
#' @param motifs Selected nucleotide contexts to be extracted.
#' @return Nucleotide context feature matrix.
#' @export

mutPredict.snv.find.motif = function(seq, motifs) {
  
  # Extract 1mer, 3mer and 5mer of all sites in roi
  seq$oneMer = substr(seq$dna, 6, 6)
  seq$threeMer = substr(seq$dna, 5, 7)
  seq$fiveMer = substr(seq$dna, 4, 8) 
  
  # Extract reverse complement for 1mer, 3mer and 5mer
  seq$onerc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq$oneMer)))
  seq$threerc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq$threeMer)))
  seq$fiverc <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(seq$fiveMer)))
  
  # Pair 1mer, 3mer and 5mer with their respective reverse complements and select A or G
  seq$one.pair = ifelse(seq$oneMer %in% c("A", "G"), seq$oneMer, seq$onerc)
  seq$three.pair = ifelse(substr(seq$threeMer, 2, 2) %in% c("A", "G"), seq$threeMer, seq$threerc)
  seq$five.pair = ifelse(substr(seq$fiveMer, 3, 3) %in% c("A", "G"), seq$fiveMer, seq$fiverc)
  
  # Extract right and left flanks for 3mer and 5mer
  seq$three.right = substr(seq$three.pair, 2, 3)
  seq$three.left = substr(seq$three.pair, 1, 2)
  seq$five.right = substr(seq$five.pair, 3, 5)
  seq$five.left = substr(seq$five.pair, 1, 3)
  
  # Assign scores to selected nucleotide context features
  seq$tM = ifelse(seq$three.pair %in% motifs$threeMer, seq$three.pair, "0")
  seq$tM = factor(seq$tM, levels = c("0", motifs$threeMer))
  seq$oM = ifelse(seq$one.pair %in% motifs$oneMer, seq$one.pair, "0")
  seq$oM = factor(seq$oM, levels = c("0", motifs$oneMer))
  seq$fM = ifelse(seq$five.pair %in% motifs$fiveMer, seq$five.pair, "0")
  seq$fM = factor(seq$fM, levels = c("0", motifs$fiveMer))
  seq$tR = ifelse(seq$three.right %in% motifs$threeRight, seq$three.right, "0")
  seq$tR = factor(seq$tR, levels = c("0", motifs$threeRight))
  seq$tL = ifelse(seq$three.left %in% motifs$threeLeft, seq$three.left, "0")
  seq$tL = factor(seq$tL, levels = c("0", motifs$threeLeft))
  seq$fR = ifelse(seq$five.right %in% motifs$fiveRight, seq$five.right, "0")
  seq$fR = factor(seq$fR, levels = c("0", motifs$fiveRight))
  seq$fL = ifelse(seq$five.left %in% motifs$fiveLeft, seq$five.left, "0")
  seq$fL = factor(seq$fL, levels = c("0", motifs$fiveLeft))
    
  seq = seq[ ,c("tM", "oM", "fM", "tR", "tL", "fR", "fL")]
  colnames(seq) = c("threeMer", "oneMer", "fiveMer", "threeRight", "threeLeft", "fiveRight", "fiveLeft")
  
  # Remove columns that are not selected
  if (length(motifs) != 1) {
    
  seq = seq[ ,names(seq) %in% names(motifs)]
  
  } else {
    
    seq = seq[names(seq) %in% names(motifs)]
    
  }
  
  return(seq)
  
}
