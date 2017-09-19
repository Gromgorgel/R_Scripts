####
# Often one needs to amplify a specific region or wants to design a primer in a narrow stretch of sequence
####
 # This function will search for the optimal primer in a short region of DNA
 # it will consider both strands and primer lengths between the specified limits
 # factors taken into account:
 # - Tm
 # - %GC
 # - G/C clamp
 # - single base repeats
 # - self complimantarity (limited)
 # The approach is purely brute-force & not intended for use on large sequence stretches !
 # the algorithm tests all possible primers of all possible lengths (within the limits)
 # therfore the number of primers it tests grows exponentially with every added base.
####
#
# This function depends on 2 other functions I've written (DNR and dnr.comp) which can be found on my github
# github.com/Gromgorgel/R_Scripts/
#
## Version History
 # V 0.0.1 = I'm just making this up as I go
 # V 0.0.2 = trying to add support for (single) degenerate bases...
 # V 0.0.3 = Tm has been made optional to allow for the construction of isothermal primers (where Tm doesn't matter much).
 #           (Tm is calcualted but does not contribute to the score, also Tm calculation has a bigger window)
 # V 0.0.4 = In an effort to improve speed of the algorithm, we switch to DNR based operations for all steps.
 #           unfortunately, that means we'll have to write some support functions to be able to get away from 
 #           string based operations.

## Known Issues
 # when counting the single base runs, the degenerate base is currently not always taken into account:
   # eg, when you have AAAM the algorithm will see it as a run of 3 bases, while in fact it is AAA A/C
   # this should be taken into account somehow. (idea: try both possibilities and take the worst score)
   # status : SOLVED
 # no support for DNAStringSet class
  # status : DEAL WITH IT (those idiots should have at least made the commands to get the sequence length the same)
 # the script is slow.
   # status : IMPROVED (I did some checks & the slowest step was BY FAR the calcualtion of the GC percentage followed by hairpin 
   #                    check (Tm calc in 3rd place) what made 'hairpin' different from the other steps was the use of 'grep'
   #                    which therefore seemed to be the bottleneck. Probably 'letterFrequency' (which previousely I used for 
   #                    GC content) also uses 'grep' or something even less efficient. I have replaced all use of 'grep' and 
   #                   'letterFrequency' by DNR based numerical approaches, so apart from a few initial steps the script
   #                    is now free of string-based operations. Elapsed time for standard inputs dropped from 15.40 to 
   #                    4.0 seconds on my system.
 # Tm window could be smaller for most usage this function will see
   # status : DEAL WITH IT (making the Tm window smaller when we're not making isothermal primers could speed up 
   #                        the function, I'll eventually get round to this)

## LEGEND
# myseq = DNAString object in which to find the primers
# limits = numeric: c(min, max) primer length
# top = numeric: number of primer candidated to return
# melt = logical: should Tm contribute to primer score
# silent =  logical: should (warning) messages be printed
# ...

optimus.primer <- function(myseq = DNAString("GAGGCAAAGCATGAAGATGATGCTGCTCTTACAGAGTTCCTTG"),
                           limits = c(18, 24), top = 10, melt = TRUE, silent = FALSE, ...){
require("DECIPHER")

# check parameters (we only continue if all boxes are checked)
################################################################################
if(class(myseq)[1] == "DNAString"){ # for the basic operations to work on the DNA sequence, it has to be in the "DNAString" format
if(length(myseq) >= limits[2]){ # seq lenght has to be at least limit
if(sum(strsplit(as.character(myseq), "")[[1]] %in% names(IUPAC_CODE_MAP)[5:15], na.rm = T) <= 1){    # do not allow more than one degenerate base

# Main function body
################################################################################
  # we make a look-up table to evaluate the parameters  calculated
  eval.tab <- matrix(c(-6, -5, -4, -3, -2, -1, 0,  1, 2, 3,                                         # Available scores
                       10,  9,  8,  7,  6,  5, 4,  3, 2, 1,                                         # bins for Tm
                       30, 25, 22, 19, 16 ,13 ,10, 6, 3, 0,                                         # bins for GC content
                       NA, NA, NA,  5,  4,  0,  1, 2, 3, NA), nrow = 4, ncol = 10, byrow = T)       # bins for GC-clamp

  # reserve space for all the primers we'll test and the scores we'll attribute
  score.table <- matrix(nrow = 2 * c(sum(length(myseq) - seq(from = limits[1], to = limits[2], by = 1))  + c(limits[2] - limits[1] + 1)), ncol = 5)
    score.table[, 1] <- rep(unlist(sapply(seq(from = limits[1], to = limits[2], by = 1), function(a){rep(a, times = length(myseq) - a + 1)})), times = 2)               # primer length
    score.table[, 2] <- rep(unlist(sapply(seq(from = limits[1], to = limits[2], by = 1), function(a){seq(from = 1, to = length(myseq) - a + 1, by = 1)})), times = 2)   # primer starting position
    score.table[, 3] <- rep(c(1, -1), each = nrow(score.table)/2)                                                                  # strand (pos / neg)
    colnames(score.table) <- c("length", "start", "strand", "Tm", "score")

  # Now that we have all the primers we want to test figured out, we can simply loop through the table and calculate the scores
  # before we start we'll message the user so that he/she has an idea of the number of primers we have to calculate
  if(!isTRUE(silent)){
    message(paste("evaluating", nrow(score.table), "possible primers", sep = " "))
    }

### let's get started
  for(i in 1:nrow(score.table)){
      ### generate the primer to be tested: depending on the strand, the positions refer to the seq or its reverse complement
        if(score.table[i, 3] == 1){
                mynuc <- myseq[score.table[i, 2]:(score.table[i, 1] + score.table[i, 2] - 1)]
        }else{
                mynuc <- reverseComplement(myseq)[score.table[i, 2]:(score.table[i, 1] + score.table[i, 2] - 1)]
        }

      ### we'll start by calculating the Tm & starting the score
        score.table[i, 4] <- seq(from = 50, to = 85, by = 0.1)[which.max(MeltDNA(DNAStringSet(mynuc), temps = seq(from = 50, to = 85, by = 0.1) , type = "derivative"))]
      if(isTRUE(melt)){
        nuc.score <- eval.tab[1, max(which(eval.tab[2,] > abs(60 - score.table[i, 4])))]
      }else{
        nuc.score <- 0
      }

      ### Now, we can do a DNR & work from that so we don't need any string based operations
        mydnr <- DNR(mynuc, type = 1)

      ### next up:GC content
        # using carefully positioned breaks, we can use histogram output to calcuate GC content
        # as a side note: degenerate bases do not count as 'full' GC points, but contribute proportional to the G/C possibility
        # eg: H ~ ACT and thus gets 1/3
        freq <- sum(hist(mydnr, breaks = c(0, 1.5, 3.5, 4.5, 13.5, 15, 23.5, 35, 123.5, 135, 235, 1500), plot = F)$counts *
                                         c(0, 1, 0, 0.5, 0, 1, 0.5, 2/3, 1/3, 2/3, 1/2))/length(mydnr) * 100
        nuc.score <- nuc.score + eval.tab[1, max(which(eval.tab[3,] > abs(50 - freq)))]

      ### add the GC-clamp score
        if(all(mydnr < 10)){
                score.update <- eval.tab[1, which(eval.tab[4,] == sum(tail(mydnr, 5) %in% c(2,3)))]
        }else{ # END no degnerate base
               # in case of a degenerate base, we'll generate all possble sequences, evaluate them all & keep only the "worst" score
                deg.pos <- which(mydnr >= 10) # which base needs replacing?
                deg.alt <- as.vector(outer(mydnr[deg.pos], 10^c(3:0), function(a,b) a%/% b %% 10)) # what are the digits of the base code
                deg.alt <- deg.alt[deg.alt != 0]# remove superluous zeros
                alt.scores <- rep(0, times = length(deg.alt))
                for(j in seq_along(deg.alt)){
                  altdnr <- c(mydnr[1:(deg.pos - 1)], deg.alt[j], mydnr[(deg.pos + 1):length(mydnr)])
                  alt.scores[j] <- eval.tab[1, which(eval.tab[4,] == sum(tail(altdnr, 5) %in% c(2,3)))]
                }# end of for loop
                score.update <- min(alt.scores) # we use a worst-case score update rather than a "mean" score
        } # END with degenerate base
        nuc.score <- nuc.score + score.update

      ### checking for single base repeats
        if(all(mydnr < 10)){
          if(any(rle(mydnr)$lengths >= 4)){ # check for runs of at least four bases
            # if there are any, we score them by their size:
            score.update <-  - 2 * (rle(mydnr)$lengths[which(rle(mydnr)$lengths > 3)] - 3.5)
            } # END of 4 base (or more) stretch
        }else{ # END no degnerate bases
            # deg.pos and deg.alt already exit. We'll reset the scores to bu sure.
            alt.scores <- rep(0, times = length(deg.alt))
            for(j in seq_along(deg.alt)){
              altdnr <- c(mydnr[1:(deg.pos - 1)], deg.alt[j], mydnr[(deg.pos + 1):length(mydnr)])
              if(any(rle(altdnr)$lengths >= 4)){ # check for runs of at least four bases
                alt.scores[j] <- - 2 * (rle(altdnr)$lengths[which(rle(altdnr)$lengths > 3)] - 3.5)
              }else{
                alt.scores[j] <- 0
              } # if-else END
            }# end of for loop
            score.update <- min(alt.scores) # we use a worst-case score update rather than a "mean" score
        } # END with degenerate base
        nuc.score <- nuc.score + score.update

      ### checking for hairpins
        # we start by creating a matrix of the positions between which we retrieve a pattern to search for in the rest of the nucleotide
        # i.e. a sliding window of size three over the first & last five bases: we'll search for the reverse complement of these, if it's present there's a chance of hairpin formation
        posmat <- matrix(c(1, 2, 3, length(mynuc) - 4, length(mynuc) - 3, length(mynuc) - 2, 3, 4, 5, length(mynuc) - 2, length(mynuc) - 1,length(mynuc)), nrow = 6, ncol = 2)
        # again we make a distinction between with/without degenerate base
        if(all(c(head(mydnr, 5), tail(mydnr, 5)) < 10)){ # since we only check the head and tail, that's where we need to check for the base
           # we start by making a list of the reverse complement of the triplets at either end of the nucleotide
           trip <- apply(posmat, 1, function(a){return(dnr.comp(mydnr[a[1]:a[2]]))})
           # check for the presence of these triplets  & attribute score depending on its position
           score.update <- sum(apply(trip, 2, function(a){
                                           sum(head(mydnr, -2) %in% a[1] & head(mydnr[-1], -1) %in% a[2] & tail(mydnr, -2) %in% a[3])}) * # returns 1 if pattern present, zero if not
                                           c(-1, -2, -3, -2, -3, -4))  # These are the scores for the self-matching. hairpin at 3' is punished harder
        }else{ # END no degnerate bases
            # deg.pos and deg.alt already exit. We'll reset the scores to bu sure.
            alt.scores <- rep(0, times = length(deg.alt))
            for(j in seq_along(deg.alt)){
              altdnr <- c(mydnr[1:(deg.pos - 1)], deg.alt[j], mydnr[(deg.pos + 1):length(mydnr)])
              trip <- apply(posmat, 1, function(a){return(dnr.comp(altdnr[a[1]:a[2]]))})
              score.update <- sum(apply(trip, 2, function(a){
                                             sum(head(altdnr, -2) %in% a[1] & head(altdnr[-1], -1) %in% a[2] & tail(altdnr, -2) %in% a[3])}) * # returns 1 if pattern present, zero if not
                                             c(-1, -2, -3, -2, -3, -4))  # These are the scores for the self-matching. hairpin at 3' is punished harder
            }# end of for loop
            score.update <- min(alt.scores) # we use a worst-case score update rather than a "mean" score
        } # END with degenerate base
        nuc.score <- nuc.score + score.update

    ## that's all the scoring for now
      score.table[i, 5] <- nuc.score
  } # END of main for loop
  ## Succesful analysis output
  trim.table <- score.table[order(score.table[, 5], decreasing = T), ][1:top, ]
  ## print the primers & info
  if(!isTRUE(silent)){
    apply(trim.table, 1, function(a){
                              if(a[3] == 1){
                              message(myseq[a[2]:(a[1] + a[2] - 1)])
                              }else{
                              message(reverseComplement(myseq)[a[2]:(a[1] + a[2] - 1)])
                              }
                              message(paste("Length:", a[1], cat("\t"), "Tm:", a[4], cat("\t"), "Strand:", a[3], cat("\t"), "Score:",  a[5], sep = " "))
                              })
    }
  return(trim.table)
################################################################################
## Error Reporting
}else{ #<=1 DEGENERATE BASES
  message("input sequence contains multiple degenerate  bases (only 1 allowed)")
  message("NA output")
}# End of degenerate check
}else{ #SEQ LENGTH
  message("input sequence shorter than upper primer length limit")
  message("NA output")
}# End of length check
}else{ #CLASS DNASTRING
  message("input sequence is not of class DNAString")
  message("NA output")
}# End of input class check

## Standard NA output
return(NA)
} # function END
