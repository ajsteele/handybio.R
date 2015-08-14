################################################################################
###  REQUIRED PACKAGES #########################################################
################################################################################

require(GenomicRanges)
require(rtracklayer)
require(ShortRead)
require(biovizBase)

# This script depends on handy.R (see https://github.com/ajsteele/handy.R), so
# make sure you source it before running this file. If you've not sourced it
# already, put it in an appropriate location and uncomment this line.
# source('handy.R')


################################################################################
###  GRANGES  ##################################################################
################################################################################

GRangesRange <- function(...) {
  # Find the range of positions encompassed by one or more GRanges objects.
  #
  # Args:
  #      ...: One or more GRanges objects.
  #
  # Returns:
  #   GRanges object spanning their range.
  grs <- list(...)
  gnm <- unique(unlist(lapply(grs, function(x){as.character(genome(x))})))
  if(length(gnm) != 1) {
    warning(paste0('The genomes differ between the GRanges ',
      'provided. GRangesRange does not work unless all sequences are from the ',
      'same genome.\ngenomes: ', paste(gnm, collapse=', ')))
  }
  sn <- unique(unlist(lapply(grs, function(x){as.character(seqnames(x))})))
  if(length(sn) != 1) {
    stop(paste0('The seqnames (eg chromosome) differ between the GRanges ',
      'provided. GRangesRange does not work unless all sequences are from the ',
      'same chromosome.\nseqnames: ', paste(sn, collapse=', ')))
  }
  min.start <- min(unlist(lapply(grs, start)))
  max.end <- max(unlist(lapply(grs, end)))
  
  GRanges(
    seqnames = sn,
    seqinfo = Seqinfo(sn, genome=gnm),
    ranges = IRanges(start = min.start,
                     end   = max.end)
  )
}

################################################################################
###  GENOMES  ##################################################################
################################################################################

liftOverPlus <- function(gr, liftover.chain.file, discard.nonunique = TRUE) {
  # Wrapper function for liftOver from rtracklayer to save some typing and deal
  # with nonunique/unliftable coordinates.
  #
  # Args:
  #      gr: GRanges object to be lifted over.
  # liftover.chain.file: string containing the name of the chain file to use.
  # discard.nonunique: logical determining whether the returned object will
  #          contain either coordinates which lift over to no locations, or
  #          multiple ones. TRUE discards, FALSE retains. (At present the zero-
  #          location coordinates are discarded regardless on flattening. Is
  #          there a better way?)
  #
  # Returns:
  #   GRanges object containing coordinates in your preferred mapping.

  # prepare for liftover! import the chain file
  liftover.chain <- import.chain(liftover.chain.file)
  # perform the liftover
  lifted.over <- liftOver(gr, liftover.chain)
  
  # if the user wants to discard those locations which don't map uniquely, then
  # do so
  if(discard.nonunique) {
    liftover.valid <- which(
      sapply(start(lifted.over), function(x){length(x) == 1}) |
      sapply(end(lifted.over),   function(x){length(x) == 1})
    )

    # return a GRanges object, with the bad ones excluded
    return(flatGrl(lifted.over[liftover.valid]))
  # if we aren't meant to discard
  } else {
    return(flatGrl(lifted.over))
    # (flatGrl discards the integer(0) list elements returned where there are no
    # matches found...is there a better way to encode their absence?)
  }
}

################################################################################
###  BIOLOGY  ##################################################################
################################################################################

revComp <- function(x) {
  # Finds the reverse complement of a DNA sequence provided as a character or
  # vector of characters.
  #
  # Args:
  #     x: A character or vector of characters containing A, C, T and G; case-
  #        insensitive.
  #
  # Returns:
  #   The reverse complement of all the elements of x.
  unlist(
    lapply(
      strsplit(x,''),# needs to be split into a vector so rev will work
      function(y) {
        paste(
          rev(
            chartr(
              # upper and lower case
              "ATGCatgc", "TACGtacg",
              y
            )
          ),
          collapse=""
        )
      }
    )
  )
}

getRestrictionFragments <- function(sequence, restriction.site) {
  # Takes a DNA sequence or sequences in a DNAStringSet and returns a GRanges
  # object containing the fragments which that sequence would be cut into.
  #
  # Args:
  #    sequence: DNAStringSet object containing a sequence or sequences.
  # restriction.site: string containing the DNA sequence bound by the
  #              restriction enzyme, including a ^ character at the point where
  #              the enzyme cuts, eg 'A^AGCTT' for HindIII (see
  #              http://en.wikipedia.org/wiki/HindIII)
  #
  # Returns:
  #   GRanges object containing restriction fragments.
  #
  # Credits:
  #   Inspired by the .getRestrictionSitesFromBSgenome function in Borbala
  #   Mifsud and Robert Sugar's GOTHiC Bioconductor package.
  #   http://www.bioconductor.org/packages/release/bioc/html/GOTHiC.html

  # find out how many base pairs from the 5' end of the restriction site the
  # enzyme cuts (subtract one because the ^ is not part of the sequence!)
  cut.pos <- strPos('^', restriction.site) - 1
  if(cut.pos == -2) {
    stop(paste0("There is no cut site (^) in your restriction.site string, '",
                restriction.site,"' (there should be exactly 1)"))
  } else if(length(cut.pos) > 1) {
    stop(paste0("There are multiple cut sites (^) in your restriction.site ",
                "string, '", restriction.site,"' (there should be exactly 1)"))
  }
  # remove ^ character from restriction site
  restriction.site <- sub('\\^', '', restriction.site)
  # find locations where restriction enzyme cuts for all chromosomes
  matches <- vmatchPattern(restriction.site, sequence)

  # return a GRanges object of the fragments
  GRanges(seqnames = rep(names(matches), sapply(matches, length) + 1),
          ranges =
        IRanges(start =
                  unlist(
                    lapply(
                      1:length(matches), # for each chromosome
                      function(i) {
                        c(1, # the first fragment starts at the chromosome start
                          # and the next ones start at the match site plus the
                          # offset to the point of actual enzyme-cutting
                          start(matches[[i]]) + cut.pos)
                      }
                    )
                  ),
                end = 
                  unlist(
                    lapply(
                      1:length(matches), # for each chromosome
                      function(i) {
                        # the last fragment ends at the end of the chromosome,
                        # and the preceding ones start at the match site plus
                        # the offset to the point of actual enzyme-cutting,
                        # plus 1 so they don't overlap...
                        c(start(matches[[i]]) + cut.pos - 1,
                          width(sequence)[i])
                      }
                    )
                  )
                )
    )
}
