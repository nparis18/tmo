Republican_vote_data.dta was generated from the raw data shipped in the
BFFPT replication package (dataverse_files.zip) by running the package's
construction chain (do-files 1_*, 1a, 1b, 1d, 1e, 1f; 1c is only needed
for appendix tables -- its permanent outputs ship pre-fixed in
Data - generated/Migration/post_lasso).
Requirements found in testing: set maxvar 32767; current ftools/reghdfe/
ivreghdfe plus require, statastates (needs an existing PERSONAL ado dir),
asgen, confirmdir, coefplot from SSC.
This is the only generated dataset needed for Table 2 of the paper and for
the tmo application in SJ/code/codeSj.do.
