### Citation:
Drokhlyansky, Smillie, van Wittenberghe, et al. (Cell, 2020)

### Description:
This file contains code for detecting heterotypic multiplets in MIRACL-Seq data. Code for other analyses found throughout the paper is available here:  
https://github.com/cssmillie/ulcerative_colitis<br/>
For other questions, please email Chris Smillie (cssmillie [AT] gmail.com).

### Requirements:
Before using this code, you must install the "dnormpar2" R package (included):  
> R CMD INSTALL dnormpar2.tar.gz

This contains a fast method for calculating the doublet statistics. 
  
### Usage:
> source('MIRACL.r')

> res = find_doublets(data, ident, sig.use)

> doublets = res[[length(res)]]$doublets

Please see code comments for more detailed usage instructions.
