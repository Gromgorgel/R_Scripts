################################################################################
# This started as a quick knock-off of the DNR function that will identify and tabulate runs of non-degenerate nucleotides
# by now it has grown so much it has deserved its own script.

## Changelog
 # 0.02 added ability to incorporate 1 flanking region (if only one degenerate base in between), set fuse to >0
 #      also added ability to select always the longest X runs of bases (set top.up to X) in case fewer pass the cutoff
 # 0.03 added support for DNAStringSet

## Known Issues
 # currently the function can only build in 1 degenerate base, ideally, the number should  be user defined
  # status : DEAL WITH IT (I have no need for more degenerate bases)
 # currently the function expects DNAString input, I have to add some lines so that DNR input works as well
  # status : DEAL WITH IT (I haven't had the time yet)

absorb.wobble <- function(myseq = DNAString("GAGGCAAAWGCATGAAGATGATGCTGCTCTTACAGSAGTTCCTTGGTGAGCAAAGCGAATCTATT"),
                  cutoff = 20, fuse = 1, top.up = 3){
require(Biostrings)
if(class(myseq)[1] == "DNAStringSet"){ # if the object class is DNAStringSet, we let the function call itself on each sequence in the set
	lapply(myseq, DNR.2)
}else{
  if(class(myseq)[1] == "DNAString"){ # for the basic operations to work on the DNA sequence, it has to be in the "DNAString" format
    # we split the DNA string into a vector of single characters
    splitt <- strsplit(as.character(myseq), "")[[1]]
    dnr <- rep(0, times = length(splitt)) # set defaults to "not ACGT" (zero)
    dnr[splitt %in% names(IUPAC_CODE_MAP)[1:4]] <- 1  # replace all ACGT positions with "1"
    # now we will use rle to locate the stretches of uninterupted 'normal' code
    dnrle <- rle(dnr)
    # Now we need to build a table that lists all stretches' start, stop and length.
    dntab <- matrix(nrow = length(dnrle$values), ncol = 5)
      dntab[, 1] <- dnrle$values
      dntab[, 2] <- dnrle$lengths
      dntab[, 3] <- c(1, head(cumsum(dnrle$lengths), -1) + 1)
      dntab[, 4] <- cumsum(dnrle$lengths)
      dntab[, 5] <- rep(0, times = length(dnrle$values))
      colnames(dntab) <- c("value", "length", "start", "end", "select")
    # For the next part we want to mark all regions that pass the cutoff & check if there are at least 3
   ## if not, we add the longest below-cutoff regions until we have at least 3 (or top.up) candidate regions
   ## quick check to see if there are indeed any stretches that pass the cut-off
    # we start by adding a mark in the select column for all above the cut-off
    dntab[dnrle$lengths >= cutoff & dnrle$values == 1, 5] <- 1
    if(sum(dntab[, 5]) < top.up){ # the number of selected regions needs topping up
        message("few non-degenerate runs longer than cutoff, topping up...")
        dntab[order(dntab[, 1], dntab[, 2] ,decreasing = T)[seq(from = 1, to = top.up, by = 1)], 5] <- 1
    } # END of top.up
    if(fuse > 0){ # extend the region to incorporate degenerate bases
        # for each selected row, we check the the region before and after (and use the longest)
        # we also check that there is only one degenerate base in between
        for(i in which(dntab[, 5] == 1)){
          # first we need a check if the region selected does not have row number 1 or 2
          # this will cause negative subscripts downstream which throws an error
          if(i <= 2){
            flanks <- if(dntab[i + 1, 2] == 1){ # check if there ar not more than 1 degenerate bases
                                                c(i + 2, dntab[i + 1, 2], dntab[i + 2, 2], 4)
                                               }else{ # if so make vector empty for downstream compatibility
                                                vector()
                                               } # END flank save
          }else if(i > length(dntab[, 1]) - 2){ # while we're at it we'll check if we don't caus a 'subscript out of bounds' error
            flanks <- if(dntab[i - 1, 2] == 1){# check if there ar not more than 1 degenerate bases
                                                c(i - 2, dntab[i - 1, 2], dntab[i - 2, 2], 3)
                                               }else{# if so make vector empty for downstream compatibility
                                                vector()
                                               } # END flank save
          }else{
            # we select the flanking regions key information (nr of bases in between, length)
            flanks <- matrix(c(i - 2, i + 2,  # the rownumbers of the next non-degenerate runs
                             dntab[c(i - 1, i + 1), 2],  # the number of degeneratge bases in between
                             dntab[c(i - 2, i + 2), 2],  # the length of the next non degenarate run
                             3, 4),                      # which column to update if this flank wins
                             nrow = 2, ncol = 4)
            # now we select the best option for fusing(if any)
            flanks <- flanks[order(flanks[, 3], decreasing = T), ] # order by length of next non degenerate run
            flanks <- flanks[ flanks[, 2] == 1, ]                  # take the top candidate after removing all those with too many degenerate bases in between
            if(length(flanks) > 4){ # check if multiple candidate remain
              flanks <- flanks[1, ]
            }
          }# END of if-else i <= 2
          # check if any flanks remain
          if(length(flanks != 0)){ #alright, we can update the start or stop of our current region
              dntab[i, flanks[4]] <- dntab[flanks[1], flanks[4]]
              dntab[i, 2] <- dntab[i, 4] - dntab[i, 3] + 1
          } # END of if != 0
        } # END of for loop
    }# END of fuse

  ## all that is left is to trim the table down to the selected rows & sort them
    dntab <- dntab[dntab[, 5] != 0, ]
    dntab <- dntab[order(dntab[, 2], decreasing = T), -c(1, 5)] # also remove none informative columns
  ## Succesful analysis output
  return(dntab)
################################################################################
  }else{ #CLASS DNASTRING
    message("input sequence is not of class DNAString")
    message("NA output")
  }# End of input class check (DNAString)
  ## Standard NA output
  return(c("length" = NA, "start" = NA, "end" = NA))
}# End of input class check (DNAStringSet)
} # function END
