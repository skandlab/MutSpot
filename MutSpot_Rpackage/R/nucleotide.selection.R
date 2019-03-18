#' Select nucleotide context through LASSO for SNVs and filter polyATs from indel mutations.
#'
#' @param sampled.sites.snv.file Sampled SNV mutated and non-mutated sites RDS file.
#' @param indel.mutations.file Indel mutations MAF file.
#' @param cutoff Frequency cutoff/threshold to determine nucleotide contexts used in prediction model, ranges from 0.5 to 1, default = 0.90.
#' @param cores Number of cores, default = 1.
#' @return A list containing frequency of nucleotide contexts selected through lasso, nucleotide contexts that passed the threshold and polyATs filtered indel mutations.
#' @export

nucleotide.selection = function(sampled.sites.snv.file, indel.mutations.file, cutoff = 0.90, cores = 1) {
  
  # If SNV mutations available, else skip this
  if (!is.null(sampled.sites.snv.file)) {
    
    print("SNV mutations available...")
    
    # Define sampled mutated and non-mutated sites
    gr = readRDS(sampled.sites.snv.file)

    # Extend to 5 bp regions
    gr2 = GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(gr)), IRanges::IRanges(IRanges::start(gr)-2, IRanges::end(gr)+2), mut = gr$mut)
    
    # Extract sequence for all sites in gr2
    extrSeq = IRanges::Views(BSgenome.Hsapiens.UCSC.hg19::Hsapiens,gr2)
    dnaSeq = GenomicRanges::as.data.frame(gr2)
    z = as(extrSeq,"DNAStringSet") 
    dnaSeq = cbind(dnaSeq, as.character(z))
    colnames(dnaSeq)[7] = "dna"
    
    # Remove sequences with N
    N.string = which(stringi::stri_detect_fixed(dnaSeq$dna, "N"))
    if (length(N.string) != 0) {
      
      N.dnaSeq = dnaSeq[N.string, ]
      dnaSeq = dnaSeq[-N.string, ]
      
    }
    
    # Extract 1mer, 3mer and 5mer
    dnaSeq$dna = as.character(dnaSeq$dna)
    dnaSeq$fiveMer = dnaSeq$dna
    dnaSeq$threeMer = substr(dnaSeq$fiveMer, 2, 4)
    dnaSeq$oneMer = substr(dnaSeq$threeMer, 2, 2) 
    
    # Pair up reverse complements
    # Reverse complement for 5mer
    dnaSeq$five.rc = Biostrings::reverse(Biostrings::chartr("ATCG", "TAGC", dnaSeq$fiveMer))
    # Reverse complement for 3mer
    dnaSeq$three.rc = Biostrings::reverse(Biostrings::chartr("ATCG", "TAGC", dnaSeq$threeMer))
    # Reverse complement for 1mer
    dnaSeq$one.rc = Biostrings::reverse(Biostrings::chartr("ATCG", "TAGC", dnaSeq$oneMer))
    
    dnaSeq$five.pair = ifelse(substr(dnaSeq$fiveMer, 3, 3) %in% c("A", "G"), dnaSeq$fiveMer, dnaSeq$five.rc)
    dnaSeq$three.pair = ifelse(substr(dnaSeq$threeMer, 2, 2) %in% c("A", "G"), dnaSeq$threeMer, dnaSeq$three.rc)
    dnaSeq$one.pair = ifelse(dnaSeq$oneMer %in% c("A", "G"), dnaSeq$oneMer, dnaSeq$one.rc)
    
    # Extract right(3') and left(5') flanks for 3mer and 5mer
    dnaSeq$three.right = substr(dnaSeq$three.pair, 2, 3)
    dnaSeq$three.left = substr(dnaSeq$three.pair, 1, 2)
    dnaSeq$five.right = substr(dnaSeq$five.pair, 3, 5)
    dnaSeq$five.left = substr(dnaSeq$five.pair, 1, 3)
    
    tab2 = aggregate(mut ~ one.pair + three.pair + three.right + three.left + five.pair + five.right + five.left, dnaSeq, FUN = function(x) { y = sum(x == 1); z = sum(x == 0); return(cbind(y, z)) }) 
    df = tab2$mut
    # Success(mutated)-1, failure(non-mutated)-0
    colnames(df) = c("Success", "Failure") 
    df = data.frame(df, tab2[ ,1:(ncol(tab2) - 1)]) 
    df = df[ ,c(3:ncol(df), 1:2)]
    
    # Covariate matrix
    one = as.factor(df$one.pair)
    if (nlevels(one) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ one - 1)
      X <- coluna
      
    }
    
    three = as.factor(df$three.pair)
    if (nlevels(three) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ three - 1)
      X <- cbind(X, coluna)
      
    }
    
    three.right = as.factor(df$three.right)
    if (nlevels(three.right) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ three.right - 1)
      X <- cbind(X, coluna)
      
    }
    
    three.left = as.factor(df$three.left)
    if (nlevels(three.left) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ three.left - 1)
      X <- cbind(X, coluna)
      
    }
    
    five = as.factor(df$five.pair)
    if (nlevels(five) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ five - 1)
      X <- cbind(X, coluna)
      
    }
    
    five.right = as.factor(df$five.right)
    if (nlevels(five.right) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ five.right - 1)
      X <- cbind(X, coluna)
      
    }
    
    five.left = as.factor(df$five.left)
    if (nlevels(five.left) > 1) {
      
      coluna <- Matrix::sparse.model.matrix( ~ five.left - 1)
      X <- cbind(X, coluna)
      
    }
    
    # Aggregated response matrix
    y_data = cbind(df$Failure, df$Success)
    colnames(y_data) = c("Failure", "Success")
    
    # Run lasso selection 100 times as a stability measure
    stabs2 = stability.sel.nucleotide(x_data = X[ ,-2], y_data = y_data, threshold = cutoff, cores = cores)
    
    feature = names(stabs2$freq.1se)
    f = stabs2$freq.1se
    freq = data.frame(feature, f)
    sel = freq
    sel = sel[order(sel$f, decreasing = TRUE), ]
    sel = as.character(sel[which(sel$f >= cutoff), "feature"])
    
  } else {
    
    freq = NULL
    sel = NULL
    
  }
  
  # If indel mutations available, else skip this
  if (!is.null(indel.mutations.file)) {
    
    print("Indel mutations available...")
    
    # Chr1-ChrX
    chrOrder <- c(paste("chr", 1:22, sep=""), "chrX")
    seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
    
    # Define indel mutations
    maf.mutations <- maf.to.granges(indel.mutations.file)
    maf.mutations = maf.mutations[as.character(GenomeInfoDb::seqnames(maf.mutations)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
    
    # Extend 9 bp on each side of the indel
    maf.mutations = maf.mutations+9
    
    # Extract sequence for all sites in maf.mutations
    extrSeq = IRanges::Views(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, maf.mutations)
    z = as(extrSeq, "DNAStringSet") 
    dnaSeq = GenomicRanges::as.data.frame(maf.mutations)
    dnaSeq$dna = as.character(z)
    
    # Remove indels with at least 8 repeated bases
    sum(grepl("AAAAAAAA",dnaSeq$dna) | grepl("TTTTTTTT",dnaSeq$dna) | grepl("CCCCCCCC",dnaSeq$dna) | grepl("GGGGGGGG",dnaSeq$dna))
    filter = which(grepl("AAAAAAAA", dnaSeq$dna) | grepl("TTTTTTTT", dnaSeq$dna) | grepl("CCCCCCCC", dnaSeq$dna) | grepl("GGGGGGGG", dnaSeq$dna))
    
    maf.mutations <- maf.to.granges(indel.mutations.file)
    if (length(filter) != 0) {
      
    df = maf.mutations[-filter]
    
    } else {
      
      df=maf.mutations
      
    }
    
    df = GenomicRanges::as.data.frame(df)
    df = df[ ,c("seqnames", "start", "end", "ral", "tal", "sid")]
    
  } else {
    
    df = NULL
    
  }
  
  return(list(freq, sel, df))
  
}





