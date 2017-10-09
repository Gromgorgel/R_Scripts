# R Scripts
A collection of assorted R scripts for Bioinformatics purposes

## the 'DNR functions.R' script
As you may have noticed, R is not really intended to handle DNA sequences. [Bioconductor](https://www.bioconductor.org/) makes many bioinormatics tools available in R. However, I found that some of the basic functions are very slow/innefficient because, well, R is not really intended to handle DNA sequences. 

Especially functions related to pattern matching/finding seem to take forever. R's built in `grep()` function, which is the most straightforward way to implement pattern matching, will slow your scripts to a crawl when used iteratively. The only way out is to **not use string-based** pattern searched but to **use number-based** pattern searches.  For shorter patterns this is rather straightforward to implement using some combination of `%in%` and logical operators. The only thing we need to do is represent the DNA as sequence of numbers rather than letters. 

enter `DNR()` or **DNA Numerical Representation** a function to translate a sequence into a numbers. There's different flavours of numerical representation to choose from (most of which I have sourced from [Mendizabal-Ruiz et al., 2017](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0173288)). 

#### integer DNR
This is the only representation to which I have made significant changes compared to what is in the above publication. The original DNR simply attributes the integers 0-3 to TCAG. Since R uses one-based indexing, I changed that to 1-4 and I changed the order of the attribution to match that of the BioStrings `IUPAC_CODE_MAP` (ACGT). Coincidentally, you can just `5 - my_dnr` to get the complement sequence. This works only if there's just standard bases because I've also added the degenerate bases (M, R, W, S, Y, K, V, H, D, B, & N) for which I have composed the numbers that represent them using the digits for the bases they encode (eg. M can be A or C so it gets the number 12). This expands the functionality of the DNR as it now functions as its own look-up table: degenerate bases can be split into their digits to yield the individual bases they encode. Also, checking for degenerate bases is as simple as `any(my_dnr > 5)` since degenerate integeres are at least 2 digits.

#### Other functions
In addition to the main `DNR()` function three more functions are contained in the script. They all offer some additional functionality when working with DNR:
- `unDNR()`: turn your number sequence back into DNA (DNAString object), currently only supports type 1 (integer DNR)
- `dnr.comp()`: turn your DNR into its complement or reverse complement (`reverse.comp = TRUE`), currently only supports type 1 (integer DNR)
- `digits()`: split a number into its digits
- `DNR.2()`: given a DNR sequence of type 1 (integer) that contains degenerate bases, this function identifies and returns the location of stretches of standard bases (eg. after doing a multiple alignment, use this to find the region where you want to design primers so that they'll amplify all sequences from the alignment).
- `dnr.explode()` : given a DNR sequence of type 1 (integer) that contains degenerate bases, this function returns a DNAStringSet of all possible non-degenerate dnr sequences that can be derived from the mother sequence.

## the 'Optimus primer.R' script
Some background: you want to amplify a specific region of some genome (because you know there are some differences there that'll let you distinguish between closely related species) and you are restricted in amplicon size because of the sequencing technology at hand, which leaves only a very narrow region of less than 40bp on either side of said region in which you can design your primers. On top of that you want to automate this so you can run it over a bunch of such regions.

Enter **optimus primer** a script that will do just that. Not an elegant script, mind you, but a purely brute-force check-all-possible-primers script. It checks Tm, %GC, G/C clamp, single base repeats, and self complimantarity (hairpins). Scores are attributed according to my personal insight & the vague-ish recommandations that float around the interwebs (I'm open to suggestions for a better scoring strategy) and out comes a list of the top performing primers. 

Users can set the desired range of primers sizes & if the primer melting temperature should contribute to the overall score. The latter has been added to allow for the design of longer isothermal PCR primers (in which melting temperature is not that important). 

'optimus primer' heavily relies on the 'DNR functions' presented above, so those have to be loaded for it to work.

## the 'run pcr.R' script
To be able to test pirmers _in silico_ a PCR simulation script is required. Some such scripts are availble throughout the different R packages, yet none managed to do exactly what I want them to do. So, building on the DNR functions I cooked up my own script. Rather than looking for perfect matches, the function takes the [PrimerMiner](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12687/abstract) approach of attributing position-dependent mismatch penalties (I have taken their scoring tables for the current implementation, so all credit belongs there). 

The main difference with PrimerMiner's _in silico_ PCR lies in the implementation: my function considers both primers at once (it looks for valid amplicons) accros the entire sequence. It also uses single DNAString objects or DNAStringSet obtjects as template rather than a sequence alignment as is the case for PrimerMiner. In this regard, 'run pcr' is more a PCR emulation, whereas the in silico PCR of 'PrimerMiner' is a primer evaluation tool.

'run pcr' heavily relies on the 'DNR functions' presented above, so those have to be loaded for it to work. In fact, once the  sequence is translated into a DNR, calculating annealing scores acros entire sequences is as straightfroward as constructing a matrix by ofsetting the sequence by one basepair per row (`nrow` = primer length), substituting each integer for its score (base integers function as their score's position in the scoring matrix), and taking the `colSum`.
