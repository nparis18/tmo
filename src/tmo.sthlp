{smcl}
{* *! version 0.9.0b1 2025-04-08}{...}
{title:Title}

{pstd}
{hi:tmo} {hline 2} Estimating standard errors via Thresholding Multiple Outcomes (TMO) method.


{title:Syntax}

{p 8 16 2}
{cmd:tmo, cmd() y() x() d()} {opth id:var()} [{it:options}]
{p_end}

{synoptset 11}{...}
{synopthdr}
{synoptline}
{synopt :{opth t:imevar(varname)}} for panel settings, the time variable must be specified in {cmd:timevar()}{p_end}
{synoptline}
{p2colreset}{...}


{marker description}{...}
{title:Description}

{pstd}
This Stata package implements the Thresholding Multiple Outcomes (TMO) method for estimating standard errors described in {help tmo##mainpaper:DellaVigna et al. (2025)}. The TMO method accounts for spatial correlation using a set of auxiliary outcomes.
 

{marker options}{...}
{title:Options}

{phang}
{opth t:imevar(varname)}} must be specified in panel settings with the time variable (e.g., year).


{marker examples}{...}
{title:Examples}

{phang}{cmd:. sysuse auto}{p_end}



{marker results}{...}
{title:Stored results}

{pstd}
{cmd:tmo} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 15 19 2: Matrices}{p_end}
{synopt:{cmd:e(tmo_se)}}TMO standard error estimate{p_end}
{p2colreset}{...}




{marker authors}{...}
{title:Author}

{pstd}
Woojin Kim{break}
NBER{break}


{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This package is in beta/testing mode. Please use at your own risk.


{marker references}{...}
{title:References}

{marker mainpaper}{...}
{phang}DellaVigna, Stefano, Guido Imbens, Woojin Kim, and David Ritzwoller. "Using Multiple Outcomes to Adjust Standard Errors for Spatial Correlation."
