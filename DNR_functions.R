# quick function to turn a sequence into a DNA Numerical Representation
 # see table 1 of Mendizabal-Ruiz et.al (2017)
 # http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0173288

## LEGEND
# type  = numeric 1 to 9. Determines numeric representation
 # 1 : integer (*)
 # 2 : real
 # 3 : EIIP
 # 4 : Atomic Number
 # 5 : Paired Numeric
 # 6 : Voss
 # 7 : Tetrahedron
 # 8 : Z- curve
 # 9 : DNA walk

 # (*) the integer DNR is the only conversion from the publication in which I have changed which
 #     value is attributed to which base. The 'new' values match the ACTG base position in the
 #     biostrings package IUPAC_CODE_MAP. Together with the fact that the numbers for degenerate
 #     bases are made up of the digits for their corresponding single bases makes the replacement
 #     of a degenerate base by its possible single bases more straightforward (ie. split intiger into its digits
 #     and use digits to look up bases in names(IUPAC_CODE_MAP) )

## Changelog
 # 0.02 added support for degenerate bases to the 'integer' DNR.
 # 0.03 added check for degenerate bases (and force "integer" if so)
 
DNR <- function(myseq = DNAString("GAGGCAAAGCATGAAGATGATGCTGCTCTTACAGAGTTCCTTG"), type = 1){
  require(Biostrings)
  if(class(myseq)[1] == "DNAString"){ # for the basic operations to work on the DNA sequence, it has to be in the "DNAString" format
    # The first five types consist of simple numerical replacements. we make a look-up table
    letterz <- names(IUPAC_CODE_MAP)
                       #     A     C     G     T       M     R     W     S     Y     K     V     H     D     B     N
    conv.tab <- matrix(c(    1,    2,    3,    4,     12,   13,   14,   23,   24,   34,  123,  124,  134,  234,  1234,  
                           -1.5,   .5,  -.5,   1.5,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,
                        0.126, 0.134, 0.0806, 0.1335, NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,
                            70,   58,   78,   66,     NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA,
                            1,   -1,   -1,     1,     NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,   NA,    NA)
                            , nrow = 5, ncol = 15, byrow = T)
    # we split the DNA string into a vector of single characters
    splitt <- strsplit(as.character(myseq), "")[[1]]
  ########################
    # run a check for degenerate bases
    if(any(splitt %in% names(IUPAC_CODE_MAP)[5:15]) & type != 1){
      message("input sequence contains degenerate bases")
      message("forcing DNR 'integer'")
      type <- 1
    }# End of input class check
  ########################
    if(type <= 5){ # putting the look-up table to work
       X1 <- rep(0, times = length(myseq))
       for(i in seq_along(letterz)){
          X1[which(splitt == letterz[i])] <- conv.tab[type, i]
          } #END of for loop
     # singular output
       splitt <- X1
      }# types 1-5 END
  ########################
    if(type == 6){ # Voss
     # prepping the components
       X1 <- rep(0, times = length(myseq))
       X2 <- rep(0, times = length(myseq))
       X3 <- rep(0, times = length(myseq))
       X4 <- rep(0, times = length(myseq))
     # saving the positions
       X1[splitt == "A"] <- 1
       X2[splitt == "G"] <- 1
       X3[splitt == "C"] <- 1
       X4[splitt == "T"] <- 1
     # singular output
       splitt <- list("X1" = X1, "X2" = X2, "X3" = X3, "X4" = X4)
    }# voss END
  ########################
    if(type == 7){ # Tetrahedron
     # prepping the components
       X1 <- rep(0, times = length(myseq))
       X2 <- rep(0, times = length(myseq))
       X3 <- rep(0, times = length(myseq))
     # saving the positions
       for(i in seq_along(splitt)){
        if(splitt[i] == "A"){
          X3[i] <- 1
        } # A end
        if(splitt[i] == "G"){
          X1[i] <- -sqrt(2)/3
          X2[i] <- -sqrt(6)/3
          X3[i] <- -1/3
        } # G end
        if(splitt[i] == "C"){
          X1[i] <- -sqrt(2)/3
          X2[i] <-  sqrt(6)/3
          X3[i] <- -1/3
        } # C end
        if(splitt[i] == "T"){
          X1[i] <- 2*sqrt(2)/3
          X3[i] <- -1/3
        } # T end
       } # for END
     # singular output
       splitt <- list("X1" = X1, "X2" = X2, "X3" = X3)
    }# tetrahedron END
  ########################
    if(type == 8){ # Z-Curve
     # prepping the components
       X1 <- rep(0, times = length(myseq))
       X2 <- rep(0, times = length(myseq))
       X3 <- rep(0, times = length(myseq))
     # initialize components
        if(splitt[1] == "A"){
          X1[1] <- -1
          X2[1] <- 1
          X3[1] <- 1
        } # A end
        if(splitt[1] == "G"){
          X1[1] <- 1
          X2[1] <- -1
          X3[1] <- -1
        } # G end
        if(splitt[1] == "C"){
          X1[1] <- -1
          X2[1] <- 1
          X3[1] <- -1
        } # C end
        if(splitt[1] == "T"){
          X1[1] <- 1
          X2[1] <- -1
          X3[1] <- 1
        } # T end
     # saving the other positions
       for(i in 2:length(splitt)){
        if(splitt[i] == "A"){
          X1[i] <- X1[i-1] - 1
          X2[i] <- X2[i-1] + 1
          X3[i] <- X3[i-1] + 1
        } # A end
        if(splitt[i] == "G"){
          X1[i] <- X1[i-1] + 1
          X2[i] <- X2[i-1] - 1
          X3[i] <- X3[i-1] - 1
        } # G end
        if(splitt[i] == "C"){
          X1[i] <- X1[i-1] - 1
          X2[i] <- X2[i-1] + 1
          X3[i] <- X3[i-1] - 1
        } # C end
        if(splitt[i] == "T"){
          X1[i] <- X1[i-1] + 1
          X2[i] <- X2[i-1] - 1
          X3[i] <- X3[i-1] + 1
        } # T end
       } # for END
     # singular output
       splitt <- list("X1" = X1, "X2" = X2, "X3" = X3)
     }# z-curve END
  ########################
    if(type == 9){ # DNA Walk
     # prepping the components
       X1 <- rep(0, times = length(myseq))
     # initialize components
        if(any(splitt[1] == c("C", "T"))){
          X1[1] <- 1
        }else{
          X1[1] <- -1
        }
     # saving the other positions
       for(i in 2:length(splitt)){
        if(any(splitt[1] == c("C", "T"))){
          X1[1] <- X1[i-1] + 1
        }else{
          X1[1] <- X1[i-1] - 1
        }# if END
       } # for END
     # singular output
       splitt <- X1
    } # dna walk END
################################################
  ## Succesful analysis output
  return(splitt)
################################################################################
  }else{ #CLASS DNASTRING
    message("input sequence is not of class DNAString")
    message("NA output")
  }# End of input class check

  ## Standard NA output
  return(NA)
} # function END


################################################################################ Another Function!
##################################################################################################
################################################################################
# quick knock-off of the above version that will identify and tabulate runs of non-degenerate nucleotides

## Changelog
 # 0.02 added ability to incorporate 1 flanking region (if only one degenerate base in between), set fuse to >0
 #      also added ability to select always the longest X runs of bases (set top.up to X)
 # 0.03 added support for DNAStingSet

DNR.2 <- function(myseq = DNAString("GAGGCAAAWGCATGAAGATGATGCTGCTCTTACAGSAGTTCCTTGGTGAGCAAAGCGAATCTATT"),
                  cutoff = 20, fuse = 1, top.up = 3){
if(class(myseq)[1] == "DNAStringSet"){ # if the object class is DNAStringSet, we let the function call itself on each sequence in the set
	lapply(myseq, DNR.2)
}else{
  require(Biostrings)
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
        dntab[order(dntab[, 1], dntab[, 2] ,decreasing = T)[seq(from = 1, to = top.up - sum(dntab[, 5]), by = 1)], 5] <- 1
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
            # now we select the best option for fusing (if any)
            flanks <- flanks[order(flanks[, 3], decreasing = T), ] # order by length of next non degenerate run
            flanks <- flanks[ flanks[, 2] == 1, ]                  # take the top candidate after removing all those with too many degenerate bases in between
            if(length(flanks) > 4){ # check if multiple candidates remain
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

################################################################################ Another Function!
##################################################################################################
################################################################################
# Since we want to abandon working with strings all together we'll write a function to make the (reverse) complement
# sequence from a type 1 DNR
dnr.comp <- function(mydnr = c(3, 1, 3, 3, 2, 1, 1, 1, 3, 2, 1, 4, 3, 1, 1, 3, 1, 4, 3, 1, 4, 3, 
			       2, 4, 3, 2, 4, 2, 4, 4, 1, 2, 1, 3, 1, 3, 4, 4, 2, 2, 4, 4, 3),
                     reverse.comp = TRUE){
  ## If there's only standard bases, we can just
  if(all(mydnr < 5)){ # there are degenerate bases, we do some more swapping
   revdnr <- 5 - mydnr
  }else{
  ## create vectors
   comple <- c(  2, 4, 3, 1, 34, 24, 14, 23, 13, 12, 234, 134, 124, 123, 1234)
   values <- c(  3, 1, 2, 4, 12, 13, 14, 23, 24, 34, 123, 124, 134, 234, 1234)
   revdnr <- rep(0, times = length(mydnr))
  # run the vectors
   for(i in seq_along(values)){
      revdnr[which(mydnr == values[i])] <- comple[i]
      } #END of for loop
  }# END of if > 5
  ## make reverse ?
  if(isTRUE(reverse.comp)){
     return(rev(revdnr))
  }else{
     return(revdnr)
  }## END of reverse
}## End of function

################################################################################ Another Function!
##################################################################################################
################################################################################
# just to be able to check the results of our above function(s) we'll write one that'll take
# a bunch of nrs and turn them into a DNA string!

unDNR <- function(mydnr = c(3, 1, 3, 3, 2, 1, 1, 1, 3, 2, 1, 4, 3, 1, 1, 3, 1, 4, 3, 1, 4, 3, 2,
                            4, 3, 2, 4, 2, 4, 4, 1, 2, 1, 3, 1, 3, 4, 4, 2, 2, 4, 4, 3)){
   # check if there aren't any 'illegal' numbers
   allowed <- c(  3, 1, 2, 4, 12, 13, 14, 23, 24, 34, 123, 124, 134, 234, 1234)
   letterz <- c("G", "A", "C", "T", "M", "R", "W", "S", "Y", "K", "V", "H", "D", "B", "N")
   if(all(mydnr %in% allowed)){
       X1 <- rep("Ni", times = length(mydnr))
       for(i in seq_along(allowed)){
          X1[which(mydnr == allowed[i])] <- letterz[i]
          } #END of for loop
################################################
  ## Succesful analysis output
   # collapsed and made into a DNAString
  return(DNAString(paste(X1, collapse = "")))
################################################################################
  }else{ #ILLEGAL NUMBERS
    message("input sequence contains illegal numbers")
    message("NA output")
  }# End of input class check

  ## Standard NA output
  return(NA)
} # function END

################################################################################ Another Function!
##################################################################################################
################################################################################
# A more elegant way to split numbers into their digits
# can be applied to vectors of numbers (will return a list)
digits <- function(x){
	if(length(x) > 1 ){
		lapply(x, digits)
	}else{
		n <- nchar(x)
		rev( x %/% 10^seq(0, length.out=n) %% 10 )
	}
}
