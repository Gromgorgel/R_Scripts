# R Scripts
A collection of assorted R scripts for Bioinformatics purposes

## the 'DNR functions.R' script
As you may have noticed, R is not really intended to handle DNA sequences. [Bioconductor](https://www.bioconductor.org/) makes many bioinormatics tools available in R. However, I found that some of the basic functions are very slow/innefficient because, well, R is not really intended to handle DNA sequences. 

Especially functions related to pattern matching/finding seem to take forever. R's built in `grep()` function, which is the most straightforward way to implement pattern matching, will slow your scripts to a crawl when used iteratively. The only way out is to **not use string-based** pattern searched but to **use number-based** pattern searches.  For shorter patterns this is rather straightforward to implement by using some combination of `%in%` and logical operators. The only thing we need to do is represent the DNA as sequence of numbers rather than letters. 

enter `DNR()` or **DNA Numerical Representation** a function to translate a sequence into a numbers. There's different flavours of numerical representation to choose from (most of which I have sourced from [Mendizabal-Ruiz et al., 2017](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0173288)). 

#### integer DNR
This is the only representation to which I have made significant changes compared to what is in the above publication. The original DNR simply attributes the integers 0-3 to TCAG. Since R uses one-based indexing, I changed that to 1-4 and I changed the order of attribution to match that of the BioStrings `IUPAC_CODE_MAP` (ACGT). Coincidentally, you can just `5 - my_dnr` to get the complement sequence. This works only if there's just standard bases because I've also added the degenerate bases (M, R, W, S, Y, K, V, H, D, B, & N) for which I have composed the numbers that represent them using the digits for the bases they encode (eg. M can be A or C so it gets the number 12). This expands the functionality of the DNR as it now functions as its own look-up table: degenerate bases can be split into their digits to yield the individual bases they encode. Also, checking for degenerate bases is as simple as `any(my_dnr > 5)` since degenerate integeres are at least 2 digits.

#### Other functions
In addition to the main `DNR()` function, four more functions are contained in the script. They all offer some additional functionality when working with DNR:
- `unDNR()`: turn your number sequence back into DNA (DNAString object), currently only supports type 1 (integer DNR)
- `dnr.comp()`: turn your DNR into its complement or reverse complement (`reverse.comp = TRUE`), currently only supports type 1 (integer DNR)
- `digits()`: split a number into its digits
- `dnr.explode()` : given a DNR sequence of type 1 (integer) that contains degenerate bases, this function returns a list of all possible non-degenerate dnr sequences that can be derived from the mother sequence.

Most of these functions have been written so they'll accept either vectors, lists of vectors, DNAStrings or DNAStringSets (depending on the fucntion). They check their input & return `NA` plus some warning message if they don't get what they want. Other than that, their use should be pretty much self-explanatory.

## the 'Absorb wobble.R' script
Given a DNA sequence that contains degenerate bases, this function identifies and returns the location of stretches of standard bases. The funtion's intended use is when performing a multiple alignments: the function can be used on the consensus sequence to find the region where one might be able to design primers so that they'll amplify all sequences from the alignment. 

The function is called 'absorb wobble' because it will fuse 2 stretches of standard bases if they are separated by exactly one degenerate base. Currenlty only a single 'absorb' is supported. In case 2 absorbs are possible (upstream and downstrea) the longer one will be picked, ties are broken at random. The function has a `cutoff` argument that sets the minimum length for a region to be considered before reporting/fusing. In addition, there is a `top.up` argument that will force the fucntion to report shorter regions if less than `top.up` regions that meet the cutoff are found.

'absorb wobble' heavily relies on the 'DNR functions' presented above, so those have to be loaded for it to work.

### Input & Output
The optimus.primer function has the following arguments:
```
absorb.wobble(myseq,  cutoff = 20, fuse = 1, top.up = 3)
```
- `myseq`: a DNAString or DNAStringSet object containing the DNA region(s) in which to search for non degenerate stretches (myseq has a standard value (not shown here) so for an example of the function output one can just `absorb.wobble()`)
- `cutoff` : numerical, non-degenerate stretches shall be at least 'cut-off' long before being considerd for absorb & reporting
- `fuse` : numerical, set to zero to disable absorbing degenerate bases, set to > 0 to enable a single absorb. In future versions this argument is intended to input the number of absorb-operations desired.
- `top.up` = numerical, if the total number of stretches that meet the 'cutoff' requirements is below 'top.up',  the next best stretches will be added until `top-up`is reached.

## the 'Optimus primer.R' script
Some background: you want to amplify a specific region of some genome (because you know there are some differences there that'll let you distinguish between closely related species) and you are restricted in amplicon size because of the sequencing technology at hand, which leaves only a very narrow region of less than 40bp on either side of said region in which you can design your primers. On top of that you want to automate this so you can run it over a bunch of such regions.

Enter **optimus primer** a script that will do just that. Not an elegant script, mind you, but a purely brute-force check-all-possible-primers script. It checks Tm, %GC, G/C clamp, single base repeats, and self complimantarity (hairpins). Scores are attributed according to my personal insight & the vague-ish recommandations that float around the interwebs (I'm open to suggestions for a better scoring strategy) and out comes a list of the top performing primers. 

Users can set the desired range of primers sizes & if the primer melting temperature should contribute to the overall score. The latter has been added to allow for the design of longer isothermal PCR primers (in which melting temperature is not that important). 

'optimus primer' heavily relies on the 'DNR functions' presented above, so those have to be loaded for it to work.

### Input & Output
The optimus.primer function has the following arguments:
```
optimus.primer(myseq, limits = c(18, 24), top = 10, melt = TRUE, silent = FALSE, ...)
```
- `myseq`: a DNAString or DNAStringSet object containing the DNA region(s) in which to design primers (myseq has a standard value (not shown here) so for an example of the function output one can just `optimus.primer()`)
- `limits`: a numerical vector of length 2, c(min, max), containing an upper and lower limit of primer length
- `top`: numerical, the number of primers to report (sorted in descending order by their score)
- `melt`: logical, should melting temperature (Tm) be taken into account when calaculating primer score? when `FALSE`Tm is calculated & reported but does not contribute to the score.
- `silent`: logical, should warning messages be reported?

The function returns a table of `top`rows and 5 columns reporting the respective length, starting position, strand (1 = forward, -1 = reverse), melting temperature, and score of each primer. When `silent = FALSE` the sequence of these primers and their main statistics are also printed on screen.

## the 'run pcr.R' script
To be able to test pirmers _in silico_ a PCR simulation script is required. Some such scripts are availble throughout the different R packages, yet none managed to do exactly what I want them to do. So, building on the DNR functions I cooked up my own script. Rather than looking for perfect matches, the function takes the [PrimerMiner](http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12687/abstract) approach of attributing position-dependent mismatch penalties (I have taken their scoring tables for the current implementation, so all credit belongs there). 

The main difference with PrimerMiner's _in silico_ PCR lies in the implementation: my function considers both primers at once accros the entire sequence (it looks for valid amplicons). It also uses single DNAString objects or DNAStringSet obtjects as template rather than a sequence alignment as is the case for PrimerMiner. In this regard, 'run pcr' is more a PCR emulation, whereas the in silico PCR of 'PrimerMiner' is a primer evaluation tool.

The current implementation of the script can handle degenerate bases (M, R, W, S, Y, K, V, H, D, B, & N) in both primer and template sequence. A 'best case scenario' is used when evaluating the match (i.e. if any of the possible bases represents a match, an exact match is assumed). However, template sequences sometimes contain long stretches of `NNNNNNNNN` due to sequencing errors or assembly problems. Such stretches of N can cause many primer matches under the 'best case scenario' evaluation, causing the algorithm to report huge numbers of amplicons, or take forever to complete. Therefore, a penalty system was incorporated that severely punishes consecutive N's (penalty increases with each additional N). The number of consecutive N that is tolerated before this procedure is called can be set by `Ntol` in the function call.

'run pcr' heavily relies on the 'DNR functions' presented above, so those have to be loaded for it to work. In fact, once the  sequence is translated into a DNR, calculating annealing scores acros entire sequences is as straightfroward as constructing a matrix by ofsetting the sequence by one basepair per row (`nrow` = primer length), substituting each integer for its score (base integers function as their score's position in the scoring matrix), and taking the `colSum`.

### Input & Output
The run.pcr function has the following arguments:
```
run.pcr(primer1, primer2, template, threshold = 100, Ntol = 5, silent = FALSE)
```
- `primer1`: a DNAString object containing the first PCR primer sequence
- `primer2`: a DNAString object containing the second PCR primer sequence
- `template`: a DNAString or DNAStringSet object containing the DNA region(s) in which the function will look for valid amplicons
- `threshold`: numerical, value above which a primer is considered not to anneal, defaults to 100 
- `Ntol`: numerical, maximum number of consecutive N in the sequence that is tolerated before penalties kick in, defaults to 5
- `silent`: logical, should warning messages be reported? In addition, if `FALSE`a progress bar will be generated when the funcion is run on a DNAStringSet

The function returns a table with one row for each primer annealing site found and with columnes `primer` (1 or 2, indicating which primer anneals at this site),  `sense` (1 or -1, fwd or reverse strand),  `position` (position of the first base of the primer of this primer binding site), `score` (annealing score, zero means perfect match), `amp_nr` (numerical identifier of the amplicon, NA if orphan primer binding site), `amp_length` (length of the amplicon, NA if orphan primer binding site).

For the function to work properly, one of the primers has to be the reverse, the other has to be the forward. The function will check both template strands for matches of the primers. All positions of the `template` are considered, mismatches stack penalties, mismatches near the 3' end receive higher penalties, as do consecutive mismatches. When `threshold` is reached the primer is considered *not* to anneal a that position. All sequence arguments have standard values (not shown here) so for an example of the function output one can just `run.pcr()` 

