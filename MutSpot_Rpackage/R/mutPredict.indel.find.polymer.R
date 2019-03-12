#' Extract polyA/T/C/G for each site in indel hotspots.
#' 
#' @param seq Dataframe containing DNA sequences of each site in the indel hotspot.
#' @return PolyA/T/C/G feature matrix.
#' @export

mutPredict.indel.find.polymer <- function(seq) {
  
  # Extract DNA sequence
  dna = seq$dna
  polymers = c("TTTTT", "AAAAA", "GGGGG", "CCCCC")
  
  # Search for polyT/A/G/C in the DNA sequences
  motifs <- sapply(polymers, function(p) {
    
    sapply(dna, function(x) {
      
      if (sum(stringr::str_detect(x, p)) >= 1)
      { 1 } else { 0 }
      
    } )
    
    } )
  colnames(motifs) = c("polyT", "polyA", "polyG", "polyC")
  
  return(motifs)
  
}
