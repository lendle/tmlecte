This R package estimates the statistical parameter `E(E(Y|A=1,B)-E(Y|A=0,B)|A=a)` using targeted maximum likelihood estimation (TMLE).  If `A` is a treatment, then this is sort of like a Treatment Effect, Conditional on one of the treatments, hence the name `tmlecte`.  Under various causal models, this statistical parameter has interpretations of things like 
 
* the affect of treatment among the treated (ATT)
* the natural direct effect (NDE), and 
* the natural direct effect among the untreated. 

The methods in this package will probably be added to the [tmle](http://cran.r-project.org/web/packages/tmle/index.html) package on CRAN in the future. 

#Installation

This package can be installed directly from github using the `devtools` package:

```R
install.packages("devtools") #if you don't have it already
library(devtools)
install_github("lendle/tmlecte")
```

Alternatively you could download the source or clone it with git, then `R CMD BUILD` and `R CMD INSTALL` it.

#Usage

Once installed, the package is loaded with `library(tmlecte)`. Documentation is available in R's built in documentation system (`?tmlecte`).  Some examples are provided in the documentation.

#References

Lendle, S. D., Subbaraman, M. S., and van der Laan, M. J. (2013). Identification and efficient estimation of the natural direct effect among the untreated. *Biometrics*, 69(2), 310-317.

Van der Laan, M. J., and Rose, S. (2011). *Targeted learning: causal inference for observational and experimental data.* Springer.
