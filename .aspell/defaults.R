## Package dictionary for R's own spell checks (utils::aspell), which CRAN's
## incoming check uses on DESCRIPTION.  inst/WORDLIST is read only by the
## 'spelling' package, not by CRAN.  The words are in nmfkc.rds:
##   saveRDS(c("FFB", "NMF", "Satoh", "Tokuda"), ".aspell/nmfkc.rds")
description <- list(encoding = "UTF-8", language = "en",
                    dictionaries = c("en_stats", "nmfkc"))
