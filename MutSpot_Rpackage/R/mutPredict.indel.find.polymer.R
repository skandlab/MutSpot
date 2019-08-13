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
  motifsA = sapply(dna, function(x) {
    
    if (sum(stringr::str_detect(x, polymers[1])) >= 1 | sum(stringr::str_detect(x, polymers[2])) >= 1) {
      
      1
      
    } else {
      
      0
      
    }
    
  } )
  
  motifsG = sapply(dna, function(x) {
    
    if (sum(stringr::str_detect(x, polymers[3])) >= 1 | sum(stringr::str_detect(x, polymers[4])) >= 1) {
      
      1
      
    } else {
      
      0
      
    }
    
  } )
  
  motifs = cbind(motifsA, motifsG)
  motifs = as.data.frame(motifs)
  colnames(motifs) = c("polyA", "polyG")
  
  return(motifs)
  
}
