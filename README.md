# R_Scripts
A collection of assorted R scripts for Bioinformatics purposes

## the 'Optimus_primer.R' script
Some background: you want to amplify a specific region of some genome (because you know there are some differences there that'll let you distinguish between closely related species) and you are restricted in amplicon size because of the sequencing technology at hand, which leaves only a very narrow region of less than 40bp on either side of said region in which you can design your primers. On top of that you want to automate this so you can run it over a bunch of such regions.

Enter **optimus primer** a script that will do just that. Not an elegant script, mind you, but a purely brute-force check-all-possible-primers script. It checks Tm, %GC, G/C clamp, single base repeats, and self complimantarity (hairpins). Scores are attributed according to my personal insight & the vague-ish recommandations that float around the interwebs (I'm open to suggestions for a better scoring strategy) and out comes a list of the top performing primers. 

Users can set the desired range of primers sizes & if the primer melting temperature should contribute to the overall score. The latter has been added to allow for the design of longer isothermal PCR primers (in which melting temperature is not that important). 

## the 'DNR_functions.R' script
As you may have noticed, R is not really intended to handle DNA sequences. [Bioconductor](https://www.bioconductor.org/), and especially its BioStrings package, makes many bioinormatics tools available in R. However, I found that some of the basic functions are very slow/innefficient because, well, R is not really intended to handle DNA sequences. 

Especially functions related to pattern matching/finding seem to take forever. R's built in `grep()` function, which is the most straightforward way to implement pattern matching, will slow your scripts to a crawl when used iteratively. The only way out is to **not use string-based** pattern searched but to **use numeric-based** pattern searches.  For shorter patterns this is rather straightforward to implement using some cobination of `%in%` and logical operators. The only thing we need to do is represent the DNA as sequence of numbers rather than letters. 

enter DNR or **DNA Numerical Representation** a function to translate a sequence into a numbers. There's different flavours of numerical representation to choose from (most of which I have sourced from [Mendizabal-Ruiz et al., 2017](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0173288)). 

### integer DNR
This is the only one representation to which I have made significant changes compared to what is in the above publication. The original DNR simply attributes the integers 0-3 to TCAG. Since R uses one-based indexing, I changed that to 1-4 and I changed the order of the attribution to match that of the BioStrings `IUPAC_CODE_MAP` (ACGT). In addition, I've added the degenerate bases (M, R, W, S, Y, K, V, H, D, B, & N) and I have compesed the numbers that represent them from the integers for the bases they encode (eg. M can be A or C so it gets the number 12). This expands the functionality of the DNR as it now functions as its own look-up table: degenerate bases can be split into their integers to yield the individual bases they encode.

### Other functions
...
