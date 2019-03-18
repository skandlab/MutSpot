#' Sample non-mutated sites and extract mutations in region of interest.
#' 
#' @param snv.mutations.file SNV mutations MAF file.
#' @param indel.mutations.file Indel mutations MAF file.
#' @param mask.regions.file Regions to mask in genome, for example, non-mappable regions/immunoglobin loci/CDS regions RDS file, default file = mask_regions.RDS.
#' @param all.sites.file All sites in whole genome RDS file, default file = all_sites.RDS.
#' @param region.of.interest Region of interest bed file, default = NULL.
#' @param ratio Sampling ratio, default = 1.
#' @param sample To sample for non-mutated sites or to use all sites in region of interest, default = TRUE.
#' @param cores Number of cores, default = 1.
#' @return A list contatining SNV mutations in region of interest, sampled mutated and non-mutated SNV sites, indel mutations in region of interest and sampled mutated and non-mutated indel sites.
#' @export

sample.sites = function(snv.mutations.file, indel.mutations.file, mask.regions.file = system.file("extdata", "mask_regions.RDS", package = "MutSpot"), all.sites.file = system.file("extdata", "all_sites.RDS", package = "MutSpot"), region.of.interest = NULL, ratio = 1, sample = TRUE, cores = 1) {

  # Chr1-ChrX
  chrOrder <- c(paste("chr", 1:22, sep=""), "chrX")
  seqi = GenomeInfoDb::seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[GenomeInfoDb::seqnames(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)[1:23]]
  
  # Define masked region i.e. CDS, immunoglobulin loci and nonmappable
  mask.regions = readRDS(mask.regions.file)
  mask.regions = mask.regions[as.character(GenomeInfoDb::seqnames(mask.regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
  
  # Define all sites in whole genome
  all.sites = readRDS(all.sites.file)
  all.sites = all.sites[as.character(GenomeInfoDb::seqnames(all.sites)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
  
  # To sample or not to sample for non-mutated sites, if the specified region is too small, may choose not to sample sites and use all sites in the specified region
  if (sample) {
    
    print("Sampling to be done...")
  
  # If snv mutations available, else skip this
  if (!is.null(snv.mutations.file)) {
    
    print("Sampling SNV sites...")
    
    # Define SNV mutations
    maf.snv.mutations <- maf.to.granges(snv.mutations.file)
    maf.snv.mutations = maf.snv.mutations[as.character(GenomeInfoDb::seqnames(maf.snv.mutations)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
    
    # If specified region, redefine SNV mutations to be in specified region
    if (!is.null(region.of.interest)) {
      
      # Define specified region
      regions = read.delim(region.of.interest, header=FALSE, sep="\t", stringsAsFactors = FALSE)
      regions = with(regions, GenomicRanges::GRanges(V1, IRanges::IRanges(V2,V3)))
      regions = GenomicRanges::reduce(regions)
      regions = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      
      ovl = IRanges::findOverlaps(maf.snv.mutations, regions)
      maf.snv.mutations = maf.snv.mutations[unique(S4Vectors::queryHits(ovl))]
      filtered.snv.mutations = maf.snv.mutations
      filtered.snv.mutations = GenomicRanges::as.data.frame(filtered.snv.mutations)
      filtered.snv.mutations = filtered.snv.mutations[ ,-c(4,5)]
      
    } else {
      
      filtered.snv.mutations = NULL
      
      }
    
    npatients.snv = length(unique(maf.snv.mutations$sid))
    maf.snv.mutations <- unique(maf.snv.mutations)

    # Remove SNV mutations in masked region
    maf.snv.mutations = subtract.regions.from.roi(maf.snv.mutations, mask.regions, cores=cores)
    
    # Target number of sites to sample, take into account of larger masked regions, and that mutated sites tend not be in masked regions
    nsites.snv = length(maf.snv.mutations)*(ratio + ratio * 0.12)
    
    # If specified region, redefine all sites to be specified region and not whole genome
    if (!is.null(region.of.interest)) {
      
      all.sites.snv = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      
      # Number of sites to sample per chromosome
      t = GenomicRanges::split(all.sites.snv, GenomeInfoDb::seqnames(all.sites.snv))
      nsites.snv.chrom = round(unlist(lapply(t, FUN=function(x) sum(as.numeric(GenomicRanges::width(x))))) / sum(unlist(lapply(t, FUN = function(x) sum(as.numeric(GenomicRanges::width(x)))))) * nsites.snv)
      seed.rand.snv = seq(1:length(t)) * 4 
      
      # Sample sites
      all.sites.snv.samples = parallel::mclapply(1:length(t), function(i) {
        
        pop = IRanges::tile(t[[i]], width = 1)
        pop = BiocGenerics::unlist(pop)
        set.seed(seed.rand.snv[i])
        pos = sample(GenomicRanges::start(pop), nsites.snv.chrom[i])
        
        if (length(pos) > 0) {
          
        gr = GenomicRanges::GRanges(unique(as.character(GenomeInfoDb::seqnames(t[[i]]))), IRanges::IRanges(pos, pos))
        return(gr)
        
        } else {
          
          return(NULL)
          
        }
        
      }, mc.cores = cores)
      
    } else {
      
      all.sites.snv = all.sites
      
      # Number of sites to sample per chromosome
      nsites.snv.chrom = round(GenomicRanges::width(all.sites.snv) / sum(as.numeric(GenomicRanges::width(all.sites.snv))) * nsites.snv)
      seed.rand.snv = seq(1:length(all.sites.snv)) * 4
      
      # Sample sites
      all.sites.snv.samples = parallel::mclapply(1:length(all.sites.snv), function(i) {
        
        set.seed(seed.rand.snv[i])
        pos = sample(GenomicRanges::start(all.sites.snv)[i]:GenomicRanges::end(all.sites.snv)[i], nsites.snv.chrom[i])
        if (length(pos) > 0) {
          
        gr = GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(all.sites.snv)[i]), IRanges::IRanges(pos, pos))
        
        } else { 
          
          return(NULL)
          
        }
        
        }, mc.cores = cores)
      
    }
    
    all.sites.snv.samples[sapply(all.sites.snv.samples, is.null)] <- NULL
    if (length(all.sites.snv.samples) == 0) {
      
      filtered.snv.mutations = NULL
      sampled.snv.sites = NULL
      
    } else {
      
    all.sites.snv.samples = suppressWarnings(do.call(getMethod(c, "GenomicRanges"), all.sites.snv.samples))
    
    # Mask selected sites that are mutated or in masked region
    mask.snv.regions = GenomicRanges::reduce(c(maf.snv.mutations, mask.regions))
    nonmut.snv.sample = subtract.regions.from.roi(all.sites.snv.samples, mask.snv.regions, cores = cores)
    nonmut.snv.sample$mut = 0
    maf.snv.mutations$mut = 1
    sampled.snv.sites = sort(c(nonmut.snv.sample, maf.snv.mutations))
    
    }
    
  } else {
    
    filtered.snv.mutations = NULL
    sampled.snv.sites = NULL
    
  }
  
  # If indel mutations available, else skip this
  if (!is.null(indel.mutations.file)) {
    
    print("Sampling indel sites...")
    
    # Define indel mutations
    maf.indel.mutations <- maf.to.granges(indel.mutations.file)
    maf.indel.mutations = maf.indel.mutations[as.character(GenomeInfoDb::seqnames(maf.indel.mutations)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
    
    # If specified region, redefine indel mutations to be in specified region
    if (!is.null(region.of.interest)) {
      
      regions = read.delim(region.of.interest ,header=FALSE, sep="\t", stringsAsFactors = FALSE)
      regions = with(regions,GenomicRanges::GRanges(V1,IRanges::IRanges(V2,V3)))
      regions = GenomicRanges::reduce(regions)
      regions = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      
      ovl = IRanges::findOverlaps(maf.indel.mutations, regions)
      maf.indel.mutations = maf.indel.mutations[unique(S4Vectors::queryHits(ovl))]
      filtered.indel.mutations = maf.indel.mutations
      filtered.indel.mutations = GenomicRanges::as.data.frame(filtered.indel.mutations)
      filtered.indel.mutations = filtered.indel.mutations[ ,-c(4,5)]
      
    } else {
      
      filtered.indel.mutations = NULL
      
      }
    
    npatients.indel = length(unique(maf.indel.mutations$sid))
    maf.indel.mutations <- unique(maf.indel.mutations)
    
    # Remove indel mutations in masked region
    maf.indel.mutations = subtract.regions.from.roi(maf.indel.mutations, mask.regions, cores = cores)
    
    # Target number of sites to sample, take into account of larger masked regions, and that mutated sites tend not be in masked regions
    nsites.indel = length(maf.indel.mutations) * (ratio + ratio * 0.12)
    
    # If specified region, redefine all sites to be specified region and not whole genome
    if (!is.null(region.of.interest)) {
      
      all.sites.indel = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      
      # Number of sites to sample per chromosome
      t = GenomicRanges::split(all.sites.indel, GenomeInfoDb::seqnames(all.sites.indel))
      nsites.indel.chrom = round(unlist(lapply(t, FUN = function(x) sum(as.numeric(GenomicRanges::width(x))))) / sum(unlist(lapply(t, FUN = function(x) sum(as.numeric(GenomicRanges::width(x)))))) * nsites.indel)
      seed.rand.indel = seq(1:length(t)) * 4    
      
      # Sample sites
      all.sites.indel.samples = parallel::mclapply(1:length(t), function(i) {
        
        pop = IRanges::tile(t[[i]], width = 1)
        pop = BiocGenerics::unlist(pop)
        set.seed(seed.rand.indel[i])
        pos = sample(GenomicRanges::start(pop), nsites.indel.chrom[i])
        if (length(pos) > 0) {
          
        gr = GenomicRanges::GRanges(unique(as.character(GenomeInfoDb::seqnames(t[[i]]))), IRanges::IRanges(pos, pos))
        return(gr)
        
        } else {
          
          return(NULL)
          
        }
        
      }, mc.cores = cores)
      
    } else {
      
      all.sites.indel = all.sites
      
      # Number of sites to sample per chromosome
      nsites.indel.chrom = round(GenomicRanges::width(all.sites.indel) / sum(as.numeric(GenomicRanges::width(all.sites.indel))) * nsites.indel)
      seed.rand.indel = seq(1:length(all.sites.indel)) * 4
      
      # Sample sites 
      all.sites.indel.samples = parallel::mclapply(1:length(all.sites.indel), function(i) {
      set.seed(seed.rand.indel[i])
      pos = sample(GenomicRanges::start(all.sites.indel)[i]:GenomicRanges::end(all.sites.indel)[i], nsites.indel.chrom[i])
      if (length(pos) > 0) {
        
      gr = GenomicRanges::GRanges(as.character(GenomeInfoDb::seqnames(all.sites.indel)[i]), IRanges::IRanges(pos, pos))
      
      } else {
        
        return(NULL)
        
      }
      
      }, mc.cores = cores)
      
    }
    
    all.sites.indel.samples[sapply(all.sites.indel.samples, is.null)] <- NULL
    if (length(all.sites.indel.samples) == 0) {
    
      filtered.indel.mutations = NULL
      sampled.indel.sites = NULL
      
    }  else {
      
    all.sites.indel.samples = suppressWarnings(do.call(getMethod(c, "GenomicRanges"), all.sites.indel.samples))
    
    # Mask selected sites that are mutated or in nonmapple regions
    mask.indel.regions = GenomicRanges::reduce(c(maf.indel.mutations, mask.regions))
    nonmut.indel.sample = subtract.regions.from.roi(all.sites.indel.samples, mask.indel.regions, cores = cores)
    nonmut.indel.sample$mut = 0
    maf.indel.mutations$mut = 1
    GenomicRanges::start(maf.indel.mutations) = GenomicRanges::start(maf.indel.mutations) + ceiling((GenomicRanges::width(maf.indel.mutations) - 1) / 2)
    GenomicRanges::end(maf.indel.mutations) = GenomicRanges::start(maf.indel.mutations)
    sampled.indel.sites = sort(c(nonmut.indel.sample, maf.indel.mutations))
    
    }
    
  } else {
    
    filtered.indel.mutations = NULL
    sampled.indel.sites = NULL
    
  }
  
  return(list(filtered.snv.mutations, sampled.snv.sites, filtered.indel.mutations, sampled.indel.sites))
    
  } else {
    
      print("No sampling...")
      
    # If SNV mutations available, else skip this
    if (!is.null(snv.mutations.file)) {
      
      print("Preparing SNV sites...")
      
      # Define SNV mutations
      maf.snv.mutations <- maf.to.granges(snv.mutations.file)
      
      # If specified region, redefine SNV mutations to be in specified region
      if (!is.null(region.of.interest)) {
        
        regions = read.delim(region.of.interest, header=FALSE, sep="\t", stringsAsFactors = FALSE)
        regions = with(regions,GenomicRanges::GRanges(V1,IRanges::IRanges(V2,V3)))
        regions = GenomicRanges::reduce(regions)
        
        ovl = IRanges::findOverlaps(maf.snv.mutations, regions)
        maf.snv.mutations = maf.snv.mutations[unique(S4Vectors::queryHits(ovl))]
        filtered.snv.mutations = maf.snv.mutations
        filtered.snv.mutations = GenomicRanges::as.data.frame(filtered.snv.mutations)
        filtered.snv.mutations = filtered.snv.mutations[ ,-c(4,5)]
        
      } else {
        
        filtered.snv.mutations = NULL
        
      }

      maf.snv.mutations = maf.snv.mutations[as.character(GenomeInfoDb::seqnames(maf.snv.mutations)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      npatients.snv = length(unique(maf.snv.mutations$sid))
      maf.snv.mutations <- unique(maf.snv.mutations)

      # Remove SNV mutations in masked region
      maf.snv.mutations = subtract.regions.from.roi(maf.snv.mutations, mask.regions, cores = cores)
      
      # If specified region, redefine all sites to be specified region and not whole genome
      if (!is.null(region.of.interest)) {
        
        all.sites.snv = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
        
      } else {
        
        all.sites.snv = all.sites
        
      }
      
      # All sites are potential non-mutated sites
      all.sites.snv = IRanges::tile(all.sites.snv, width = 1)
      all.sites.snv = BiocGenerics::unlist(all.sites.snv)

      # Mask selected sites that are mutated or in masked region
      mask.snv.regions = GenomicRanges::reduce(c(maf.snv.mutations, mask.regions))
      nonmut.snv.sample = subtract.regions.from.roi(all.sites.snv, mask.snv.regions, cores = cores)
      nonmut.snv.sample$mut = 0
      maf.snv.mutations$mut = 1
      sampled.snv.sites = sort(c(nonmut.snv.sample, maf.snv.mutations))
      
    } else {
      
      filtered.snv.mutations = NULL
      sampled.snv.sites = NULL
      
    }
    
    # If indel mutations available, else skip this
    if (!is.null(indel.mutations.file)) {
      
      print("Preparing indel sites...")
      
      # Define indel mutations
      maf.indel.mutations <- maf.to.granges(indel.mutations.file)
      
      # If specified region, redefine indel mutations to be in specified region
      if (!is.null(region.of.interest)) {
        
        regions = read.delim(region.of.interest,header=FALSE, sep="\t", stringsAsFactors = FALSE)
        regions = with(regions,GenomicRanges::GRanges(V1,IRanges::IRanges(V2,V3)))
        regions = GenomicRanges::reduce(regions)
        
        ovl = IRanges::findOverlaps(maf.indel.mutations, regions)
        maf.indel.mutations = maf.indel.mutations[unique(S4Vectors::queryHits(ovl))]
        filtered.indel.mutations = maf.indel.mutations
        filtered.indel.mutations = GenomicRanges::as.data.frame(filtered.indel.mutations)
        filtered.indel.mutations = filtered.indel.mutations[ ,-c(4,5)]
        
      } else {
        
        filtered.indel.mutations = NULL
        
      }
      
      maf.indel.mutations = maf.indel.mutations[as.character(GenomeInfoDb::seqnames(maf.indel.mutations)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
      npatients.indel = length(unique(maf.indel.mutations$sid))
      maf.indel.mutations <- unique(maf.indel.mutations)
      
      # Remove indel mutations in masked region
      maf.indel.mutations = subtract.regions.from.roi(maf.indel.mutations, mask.regions, cores = cores)
      
      # If specified region, redefine all sites to be specified region and not whole genome
      if (!is.null(region.of.interest)) {
        
        all.sites.indel = regions[as.character(GenomeInfoDb::seqnames(regions)) %in% as.character(GenomeInfoDb::seqnames(seqi))]
        
      } else {
        
        all.sites.indel = all.sites
        
      }
      
      # All sites are potential non-mutated sites
      all.sites.indel = IRanges::tile(all.sites.indel, width = 1)
      all.sites.indel = BiocGenerics::unlist(all.sites.indel)
        
        
      # Mask selected sites that are mutated or in masked region
      mask.indel.regions = GenomicRanges::reduce(c(maf.indel.mutations, mask.regions))
      nonmut.indel.sample = subtract.regions.from.roi(all.sites.indel, mask.indel.regions, cores = cores)
      nonmut.indel.sample$mut = 0
      maf.indel.mutations$mut = 1
      GenomicRanges::start(maf.indel.mutations) = GenomicRanges::start(maf.indel.mutations) + ceiling((GenomicRanges::width(maf.indel.mutations) - 1) / 2)
      GenomicRanges::end(maf.indel.mutations) = GenomicRanges::start(maf.indel.mutations)
      sampled.indel.sites = sort(c(nonmut.indel.sample, maf.indel.mutations))
      
    } else {
      
      filtered.indel.mutations = NULL
      sampled.indel.sites = NULL
      
    }
    
    return(list(filtered.snv.mutations, sampled.snv.sites, filtered.indel.mutations, sampled.indel.sites))
    
    }
}
