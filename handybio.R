################################################################################
###  REQUIRED PACKAGES #########################################################
################################################################################

require(GenomicRanges)
source('handy.R')


################################################################################
###  GRANGES  ##################################################################
################################################################################

dataFrame2GRanges <- function(df, excludeCols = NULL, keepMeta = TRUE,
                              keepStrand = ("strand" %in% names(df))) {
  # Samples rows from a data frame.
  #
  # Args:
  #      df: A data frame containing data to be converted to a GRanges object.
  #          Certain columns are mandatory: at least two of start, end, width;
  #          one of chr, seqnames; strand if keepStrand is TRUE
  # excludeCols: A vector of column names or numbers to exclude from the GRanges
  #          object generated.
  # keepMeta: Boolean, whether or not to keep additional columns in the
  #          data frame as metadata in the generated GRanges object.
  # keepStrand: Boolean, whether or not to keep strand information. Defaults to
  #          TRUE if a column called strand is present in the data frame.
  #
  # Returns:
  #   GRanges object generated from the data frame.
  #
  # Credits:
  #   Inspired by Kasper Daniel Hansen's function data.frame2GRanges at
  #   https://stat.ethz.ch/pipermail/bioconductor/2011-November/042333.html
  
  # if df isn't a data frame
  if(class(df) != "data.frame") {
    stop(paste0(deparse(substitute(df)), ' is not a data frame.'))
  }
  # if there aren't two of start, end or width, throw an error
  if (sum(c("start", "end", "width") %in% names(df)) < 2) {
    if (sum(c("start", "end", "width") %in% names(df)) < 2) {
      error.info <- 'you have none!'
      } else {
        error.info <- paste('you have only',
          c("start", "end", "width")[c("start", "end", "width") %in% names(df)])
      }
    stop(paste0('You need at least two of: start, end, width'))
  }
  # if the chromosome names aren't defined somewhere, throw an error
  if (!any(c("chr", "seqnames") %in% names(df))) {
      stop('You need to include chromosome name as a column called chr or ',
           'seqnames')
  # if they are, make sure naming is consistent with the rest of the script
  } else if ("seqnames" %in% names(df)) {
      names(df)[names(df) == "seqnames"] <- "chr"
  }
  # if we're keeping strand data, make sure it's in the +/- format preferred
  # by GRanges objects
  if(keepStrand && "strand" %in% names(df)) {
    if(is.numeric(df$strand)) {
      strand <- ifelse(df$strand == 1, "+", "*")
      strand[df$strand == -1] <- "-"
      df$strand <- strand
    }
  }
  # build the GRange
  gr <- GRanges(seqnames = df$chr,
                ranges = IRanges(start = df$start,
                                 end   = df$end,
                                 width = df$width),
                strand = (if(keepStrand) df$strand else "*")
                )
  # if we're keeping metadata, append those columns not already used to the
  # GRanges object as metadata
  if(keepMeta) {
    mcols(gr) <- df[ , 
     !(names(df) %in% c("chr", "start", "end", "width", "strand", excludeCols)),
     drop = FALSE # return a data frame even if only one column--preserves name
                   ]
  }
  # finally, make the rownames consistent
  names(gr) <- rownames(df)
  # and return the GRanges object
  gr
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