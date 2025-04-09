# `TMO`:  Thresholding Multiple Outcomes (TMO)
Stata package to implement TMO, which uses multiple outcomes to adjust standard errors for spatial correlation.

PACKAGE IS CURRENTLY FOR TESTING ONLY.

# Install
To install the TMO package from this repository, please copy and run the following lines in Stata:
<pre> * Remove existing program if installed previously
  cap ado uninstall tmo
  
  * Install most recent version
  net install tmo, from("https://raw.githubusercontent.com/wjnkim/tmo/master/src")
</pre>

# References
DellaVigna, Stefano, Guido Imbens, Woojin Kim, and David Ritzwoller. (2025). "Using Multiple Outcomes to Adjust Standard Errors for Spatial Correlation."
