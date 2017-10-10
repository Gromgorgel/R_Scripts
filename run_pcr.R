# This in silico PCR function follows the scoring strategy from PrimerMiner as well as their score & penalty tables
# "but why don't you just use PrimerMiner??" I hear you ask
# well, this functions is more flexible vis-a-vis my current workflow and data.
# Also, it considers both primers at once, it tests all possibilities for annealing in a sequence & it works on DNAStrings 
# rather than on multiple sequence alignments.
# so there's that.

## Changelog
 # 0.02 added support for DNAStringSet

## Known Issues
 # currently the function can only deal with 1 degenerate base per primer, ideally, the number shouldn't matter
  # status : DEAL WITH IT (I have no need for more degenerate bases)
 # currently the function cannot deal with degenerates in the template
  # status : DEAL WITH IT (you can present the different possibilities for the template as a DNAStringSet,
  #                        also see dnr.explode() with argument 'undnr = TRUE' in DNR_functions)

# for debugging purposes:
#  templateSet <- DNAStringSet(c("seq1" = "TTTTACTGCCAACCAAGGATGTCAGATGATATCCTTGGTTGGCACTAGGATATGGATCTTCATGGAAAGGTCGTGTGGGAGACACTTTTGGATGCTGAA",
#                                "seq2" = "TTTTACTGCCAACCAAGGATGTCAGATGATATCCTTGGTTGGCACTTATGGATCTTCATGGAAAGGTCGTGTGGGAGACACTTTTGGATGCTGAA"))

run.pcr <- function(primer1 = DNAString("TACTGCCAACCAAGGATRTCA"),
                    primer2 = DNAString("GCATCCAAAAGTGTCTCCCA"),
                    template = DNAString("TTTTACTGCCAACCAAGGATGTCAGATGATATCCTTGGTTGGCACTAGGATATGGATCTTCATGGAAAGGTCGTGTGGGAGACACTTTTGGATGCTGAA"),
                    threshold = 100 ){ # score below which amplification is considered impossible

# check parameters (we only continue if all boxes are checked)
################################################################################
require(Biostrings)
## check for StringSet
if(class(template)[1] == "DNAStringSet"){ # if the object class is DNAStringSet, we let the function call itself on each sequence in the set
  # it's slightly more complicated than for our other functions since we have to pass each primer to the different sequences in the set
  namnam <- names(template)
  if(is.null(namnam)){ # no names! let's make some
    namnam <- apply(matrix(c(rep("seq", times =  length(templateSet)), seq_along(templateSet)),
              nrow = length(templateSet), ncol = 2), 1, function(a){paste(a[1], a[2], sep = "")})
  } #END of names
  amp.Set.table <- list()
  for(x in seq_along(template)){
      amp.Set.table[[namnam[x]]] <- run_pcr(primer1, primer2, template[[x]])
  } # END of x for

}else{
## Check for DNAString
if(class(template)[1] == "DNAString"){ # for the basic operations to work on the DNA sequence, it has to be in the "DNAString" format

## check for Degenerate bases
  # first we transform the letter sequence into a number sequence
  dnr.pr1 <- DNR(primer1,  type = 1)
  dnr.pr2 <- DNR(primer2,  type = 1)
  dnr.tpl <- DNR(template, type = 1)
  # we allow 1 degenerate base per primer
if(sum(dnr.pr1 > 4) < 2 & sum(dnr.pr2 > 4) < 2){
  # We currently do NOT allow degenerate bases in the template
  # now we add a check for the number of degenerate bases:
if(sum(dnr.tpl > 4) == 0){

# Main function body
################################################################################
## Scores & penalties are taken from: Elbrecht & Leese, Methods in Ecology and Evolution 2017, 8, 622-626

  ## score conversion matrix 
   # by having the row & col numbers correspond to the base DNR integers we can select the right columns like: score.conv[dnr.pr2,]
  score.conv <-  matrix(c(  0,   1, 0.5,   1,  #  A = 1
                            2,	 0,   2,  0.5,  #  C = 2
                          0.5,   2,   0,   1,  #  G = 3
                            2, 0.5,   2,	  0), #  T = 4
                        nrow = 4, ncol = 4, byrow = T)
  pos.penalty <- c(242.4, 202.8, 169.8, 142.4, 119.5, 100.4, 84.5 , 71.2, 60.2, 51, 43.3, 36.9,
                   31.6, 27.2, 23.5, 20.4, 17.8, 15.7, 13.9, 12.4, 11.2, 10.2, 9.3, 8.6, 8,
                   7.5, 7.1, 6.7, 6.4, 6.2)

  # Looper
  ##############################################################################
  ## we will do a complete analysis: match both primers in both directions,
   # to do so, we need to create 2 loops: one over the primers (p), one over the directions (d)
   # we also create a loop to accomodate any degenerate bases in the primer
  for(d in 1:2){ ## direction loop

    for(p in 1:2){  ## Primer loop
      primer <- get(paste("dnr.pr", p, sep = ""))
       ## set direction:
        if(d == 2){
           primer <- dnr.comp(primer)
        }# END if d

      # this is also where we create the binding position matrix so we can add to it (in case of degenerates)
      binder <- matrix(c(NA, NA), nrow = 2, ncol = 1)

      ## check for degenerate bases
      if(any(primer > 4)){
        wub.pos <- which(primer > 4)
        wub <- digits(primer[wub.pos])
      }else{ # no wobbles
        wub <- 0
      }# END of

      for(w in wub){ ## wobble loop
        if(w != 0){  # if there are wobbles, we replace them by their digits
          primer[wub.pos] <- w
        } # END if w

        # Actual primer matching
        ########################################################################
        ## we construct the matching matrix
         # the matching matrix consists of 'nrow' repeats of the sequence, each next repeat has one base of the front less than the previous row
         # in essence, each column represents an annealing position of the primer (ie. matchbox[, 1] are ther first 'nrow' bases etc)
         # if we now score each template-to-primer match/mismatch the colsums will tell us how well the primer anneals at that position
         # due to the nature of the scores (match = zero, mismatches stack penalties) if colsum == zero, we have a perfect match.
        matchbox <- matrix(c(rep(dnr.tpl, times = length(primer)), rep(NA, times = length(primer))),
                           nrow = length(primer), ncol = length(dnr.tpl) + 1, byrow = TRUE)
        # cut-off the last 'nrow' columns as these are not informative (primer overhang)
        matchbox <- matchbox[ , head(1:dim(matchbox)[2], -length(primer))]

        ## create empty score matrix
        scorebox <- matrix(nrow = length(primer), ncol = dim(matchbox)[2])
        # as the dnr functions as its own look-up we can just use the base integers to retrieve their match score
        for(i in seq_along(primer)){ # loop primer
           for(j in 1:dim(matchbox)[2]){ # loop sequence
              scorebox[i, j] <- score.conv[primer[i], matchbox[i, j]]
           } # END of j loop
        } # END of i loop

        ## Next we multiply each score with a penalty depending on its position in the primer
        scorebox <- sweep(scorebox, MARGIN = 1, pos.penalty[1:length(primer)], '*')

        ## to complete the score, all we need to do is increase the penalty scores for adjecent mismatches
        for(i in 1:dim(scorebox)[2]){ # loop sequence
           for(j in 1:(length(primer) - 1)){ # loop primer
             scorebox[j:(j+1), i] <- scorebox[j:(j+1), i] * sum(scorebox[j:(j+1), i] != 0)
             # the above just multiplies the penalties times 2 if both are a mismatch, times 1 if one is a mismatch
             # and times zero if both are a match (but then they were already zero to begin with)
             # note that when a mismatch is flanked by two other mismatches it thus get multiplied by 4
             # from what I can tell this is also the case in the PrimerMiner code.
           } # END of j loop
        } # END of i loop

        bind <- rbind(c(1:dim(scorebox)[2]), colSums(scorebox))[, colSums(scorebox) <= threshold]
        # make sure we have a matrix (if only one column R simplifies to a vector)
        if(length(bind) == 2){
           bind <- matrix(bind, nrow = 2, ncol = 1)
        }else if(length(bind) == 0){
           bind <- matrix(c(NA, NA), nrow = 2, ncol = 1)
        }# END if else
        # we add the bininding sites to the binder vector, except if the binding site is already
        # there (from another wobble integer), in that case we update the annealing score (if it's lower & thus more likely)
        if(any(bind[1, !is.na(bind[1,])] %in% binder[1, !is.na(binder[1,])])){ # some binding positions are already present (excluding NAs)
          dupes <- unique(c(bind[1, ], binder[1, ])[duplicated(c(bind[1, ], binder[1, ]))]) # duplicated values
          dupes <- dupes[!is.na(dupes)] # remove NA values from duplicates
          # we start by adding the non-duplicated ones to binder
          binder <- cbind(binder, bind[,!bind[1, ] %in% dupes])
          for(v in dupes){ # loop through duplicated values
            binder[2, binder[1, ] == v] <- min(c(binder[2, binder[1, ] == v], bind[2, bind[1, ] == v]), na.rm = TRUE)
          } # END for v
        }else{ # no duplicated values
          binder <- cbind(binder, bind)
        } # END any duplicated
      } # END for w

    ## we now save the results for the current primer & direction comibination (before going to the next primer)
    varname <-  paste("pr", p, "_", d, sep = "")
    assign(varname, binder)
    } # END for p
  } # ENd for d

  # somewhere we have to cbind the "bind" matrices so that the wobbles are joined
  # by creating an empty vector at the start this should work out when there are no wobbles as well
  amp.table <- matrix(c("primer" = NA, "sense" = NA, "position" = NA, "score" = NA, "amp_nr" = NA, "amp_length" = NA),
                      nrow = 1, ncol = 6)
  amp.counter <- 0
  ## checking for amplicons and building the annealing table
  if(any(!is.na(pr1_1[1, ]))){ # any matches for pr1 in fwd?
    sites1 <-  pr1_1[1, !is.na(pr1_1[1, ])] # list of matches for pr1
    for(site1 in sites1){
      if(any(pr2_2[1, !is.na(pr2_2[1, ])] > site1)){ # is there any non-NA binding site for pr2 in rev?
         sites2 <-  pr2_2[1, !is.na(pr2_2[1, ])] # list of matches for pr2
         for(site2 in sites2){
           if(site2 > site1){ # we have amplicon!
            amp.counter <- amp.counter + 1 # step counter
            amp <- matrix(c(1,  1, site1, pr1_1[2, which(pr1_1[1, ] == site1)], amp.counter, site2 - site1 + length(dnr.pr2),
                            2, -1, site2, pr2_2[2, which(pr2_2[1, ] == site2)], amp.counter, site2 - site1 + length(dnr.pr2)),
                            nrow = 2, ncol = 6, byrow = T)
            colnames(amp) <- c("primer", "sense", "position", "score", "amp_nr", "amp_length")
            # AND we add it to the amp.table
            amp.table <- rbind(amp.table, amp)
           } # END site2 > site1
         }# END for site2
      }# END if pr2_2
    }# END for site1
  }# END if pr1_1

  # we now do the same for amp2 F & amp1 R
  if(any(!is.na(pr2_1[1, ]))){ # any matches for pr2 in fwd?
    sites1 <-  pr2_1[1, !is.na(pr2_1[1, ])] # list of matches for pr2
    for(site1 in sites1){
      if(any(pr1_2[1, !is.na(pr1_2[1, ])] > site1)){ # is there any non-NA binding site for pr1 in rev?
         sites2 <-  pr1_2[1, !is.na(pr1_2[1, ])] # list of matches for pr2
         for(site2 in sites2){
           if(site2 > site1){ # we have amplicon!
            amp.counter <- amp.counter + 1 # step counter
            amp <- matrix(c(1,  1, site1, pr2_1[2, pr2_1[1, ] == site1], amp.counter, site2 - site1 + length(dnr.pr2),
                            2, -1, site2, pr1_2[2, pr1_2[1, ] == site2], amp.counter, site2 - site1 + length(dnr.pr2)),
                            nrow = 2, ncol = 6, byrow = T)
            colnames(amp) <- c("primer", "sense", "position", "score", "amp_nr", "amp_length")
            # AND we add it to the amp.table
            amp.table <- rbind(amp.table, amp)
           } # END site2 > site1
         }# END for site2
      }# END if pr2_2
    }# END for site1
  }# END if pr1_1

  # we also add the other annealing sites to the table (if any)
  for(d in 1:2){ ## direction loop
    for(p in 1:2){  ## Primer loop
    pos <- get(paste("pr", p, "_", d, sep = ""))
    to.add <- which(!pos[1, ] %in% amp.table[, 3])
    nr <- length(pos[1, to.add])
    sns <- if(d == 1){ 1 }else{ -1 }
    amp.table <- rbind(amp.table, matrix(c(rep(p, times = nr), rep(sns, times = nr), pos[1, to.add], pos[2, to.add], 
                                           rep(NA, times = nr), rep(NA, times = nr)), ncol = 6, nrow = nr, byrow = F))
    } # END for p
  } # ENd for d

 # in a last step, we check if the amp.table has grown (if so we cut off the first line which were placeholder NAs)
 if(length(amp.table) > 5){
   amp.table <- amp.table[-1, ]
   } # END length amp.table

 map <- setNames(c("fwd", "rev"), c(1,-1))
 rownames(amp.table) <- map[as.character(amp.table[, 2])]
 
 
  ## Succesful analysis output
  return(amp.table)

################################################################################
## Error Reporting
  }else{ # WOBBLY TEMPLATE
    message("template sequence contains degenerate bases")
    message("NA output")
  }# End of input check for degenerate bases (template)

  }else{ # WOBBLY PRIMERS
    message("primer sequences contain too many degenerate bases")
    message("NA output")
  }# End of input check for degenerate bases (primers)

  }else{ #CLASS DNASTRING
    message("template sequence is not of class DNAString")
    message("NA output")
  }# End of input class check (DNAString)

  ## Standard NA output
    return(c("primer" = NA, "sense" = NA, "position" = NA, "score" = NA, "amp_nr" = NA, "amp_length" = NA))
}# End of input class check (DNAStringSet)
  return(amp.Set.table)
} # function END
