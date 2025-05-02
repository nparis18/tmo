# `TMO`:  Thresholding Multiple Outcomes (TMO)
Stata package to implement TMO, which uses multiple outcomes to adjust standard errors for spatial correlation.
[Stata documentation](https://github.com/wjnkim/tmo/blob/main/tmo.pdf)

PACKAGE IS CURRENTLY FOR TESTING ONLY.

# Installation
To install the TMO package from this repository, please copy and run the following lines in Stata:
<pre> * Remove existing program if installed previously
  cap ado uninstall tmo
  
  * Install most recent version
  net install tmo, from("https://raw.githubusercontent.com/wjnkim/tmo/master/src")
</pre>

# References
DellaVigna, Stefano, Guido Imbens, Woojin Kim, and David Ritzwoller. (2025). "Using Multiple Outcomes to Adjust Standard Errors for Spatial Correlation." [https://arxiv.org/pdf/2504.13295](https://arxiv.org/pdf/2504.13295)
