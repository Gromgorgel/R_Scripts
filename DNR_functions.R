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
    conv.tab <- matrix(c(    1,    2,    3,    4,     12,   13,   14,   23,   24,   34,  123,  124,  134,  234,  1234,  # note that the digits of the integers for the degenerate bases indicate which bases they represent
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
