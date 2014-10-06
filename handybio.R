################################################################################
###  REQUIRED PACKAGES #########################################################
################################################################################

require(GenomicRanges)


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
  
  # if df isn't a data frame,
  if(class(df) != "data.frame") {
    stop(paste0(deparse(substitute(df)), ' is not a data frame.'))
  }
  # if there aren't two of start, end or width, throw an error
  if (sum(c("start", "end", "width") %in% names(df)) < 2) {
      stop('You need at least two of: start, end, width')
  }
  # if the chromosome names aren't defined somewhere, throw an error
  if (!any(c("chr", "seqnames") %in% names(df))) {
      stop('You need to include chromosome name as chr of seqnames')
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
                strand = ifelse(keepStrand, df$strand, "*")
                )
  # if we're keeping metadata, append those columns not already used to the
  # GRanges object as metadata
  if(keepMeta) {
    elementMetadata(gr) <- df[ , 
        !(names(df) %in% c("chr", "start", "end", "width", "strand",
                           excludeCols)
        )
        ]
  }
  # finally, make the rownames consistent
  names(gr) <- rownames(df)
  # and return the GRanges object
  gr
}