{smcl}
{* *! version 0.9.0b1 2025-04-08}{...}
{title:Title}

{pstd}
{hi:tmo} {hline 2} Estimating standard errors via the Thresholding Multiple Outcomes (TMO) method.

{title:Syntax}

{p 8 16 2}
{cmd:tmo}, {cmd:cmd}({it:cmdline}) {cmd:x}({it:{help varname}}) {cmd:ylist}({it:{help varlist}}) {cmd:{opt i:dvar}}({it:{help varname}}) [{help tmo##options_table:options}]
{p_end}

{marker options_table}{...}
{synoptset 22 tabbed}{...}
{synopthdr}
{synoptline}
{synopt :{opt cmd(cmdline)}} {it:cmdline} is the command that produces the 
regression of interest{p_end}
{synopt : } - {cmd:tmo} currently supports regressions using {help regress}, 
{help reghdfe}, {help ivreg2}, or {help ivreghdfe}{p_end}

{synopt :{opth x(varname)}} regressor of interest in {it:cmdline} for which to 
estimate TMO standard errors{p_end}
{synopt : } - {cmd:tmo} estimates the standard error for only this declared 
independent variable {p_end}

{synopt :{opth ylist(varlist)}} list of auxiliary outcomes to use in {cmd:tmo} 
{p_end}

{synopt :{opth i:dvar(varname)}} location identifier variable; must be unique 
(within each {it:t} for panel case) {p_end}

{syntab:Panel setting}
{synopt :{opth t:imevar(varname)}} time identifier variable; must be declared 
for panel case{p_end}

{syntab:Optional}
{synopt :{opt miss:limit(#)}} limit for proportion of observations allowed to be 
missing for auxiliary outcomes {p_end}
{synopt : } - auxiliary outcomes missing more than {cmd:misslimit} are not used 
{p_end}
{synopt : } - {cmd:misslimit} must be in [0,1]; default is 0.1 {p_end}

{synopt :{opt file:suffix(str)}} folder path and base filename for saving 
figures and results {p_end}
{synopt : } - required for {cmd:plot} or {cmd:save} options below{p_end}

{synopt :{opt savedyad}} save Stata data file with correlation and contribution 
to standard error for each location pair {p_end}

{synopt :{opt plotq}} save plot of optimal threshold estimator {p_end}

{synopt :{opt plothist}} save plot for histogram of correlations between 
locations {p_end}
{synopt :{opt plothistnbins(#)}} number of bins for histogram of correlations 
(default 10000) {p_end}

{synopt :{opt plotse}} save plot for standard error estimates across thresholds 
{p_end}
{synopt :{opt saveplotseest}} save Stata data file with standard error estimates 
across thresholds {p_end}

{synopt :{opt saveest}} save results in {cmd:r()} to Stata data file {p_end}
{synoptline}
{p2colreset}{...}

{marker description}{...}
{title:Description}

{pstd}
This Stata package implements the Thresholding Multiple Outcomes (TMO) method 
for estimating standard errors in {help tmo##paper:DellaVigna et al. (2025)}. 
The TMO method accounts for spatial correlation between locations by using a 
set of auxiliary outcomes to estimate the correlation of errors between 
locations. TMO allows pairs of locations with correlations above the optimal 
threshold to have correlated error terms in the standard error estimate.

{pstd}
To use the {cmd:tmo} package, enter the Stata command that produces the 
regression of interest in the {cmd:cmd()} option. {cmd:tmo} calculates the TMO
standard error for the independent variable specified in {cmd:x()}.

{marker examples}{...}
{title:Examples}

{pstd}Load US county dataset{p_end}
{phang2}{cmd:. use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_differences.dta", clear}{p_end}

{pstd}Define list of auxiliary outcomes{p_end}
{phang2}{cmd:. qui ds fips stfips PIN_persincpc_d EDU_college_d, not}{p_end}
{phang2}{cmd:. local ylist `r(varlist)'}{p_end}

{pstd}Regression of interest: change in per capita income on change in college 
educated with state fixed effects{p_end}
{phang2}{cmd:. reg PIN_persincpc_d EDU_college_d i.stfips, r}{p_end}

{pstd}Run TMO{p_end}
{phang2}{cmd:. tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips)}{p_end}

{pstd}Run TMO and save figures{p_end}
{phang2}{cmd:. tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips) file(./example) plotq plothist plotse}{p_end}
{hline}

{pstd}Panel example{p_end}
{phang2}{cmd:. use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_panel.dta", clear}{p_end}
{phang2}{cmd:. qui ds fips stfips EMN_farm EDU_publicenroll year, not}{p_end}
{phang2}{cmd:. local ylist `r(varlist)'}{p_end}
{phang2}{cmd:. tmo, cmd(reghdfe EMN_farm EDU_publicenroll i.year, absorb(stfips) cluster(fips)) x(EDU_publicenroll) ylist(`ylist') i(fips) t(year) file(./example_panel) plotq plothist}{p_end}

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:tmo} stores the following in {cmd:r()}:

{synoptset 22 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:r(beta)}} coefficient on {cmd:x()}{p_end}
{synopt:{cmd:r(tmo_se)}} TMO standard error estimate{p_end}
{synopt:{cmd:r(orig_se)}} original standard error from {cmd:cmd()}{p_end}
{synopt:{cmd:r(lb)}} lower bound of 95% confidence interval{p_end}
{synopt:{cmd:r(ub)}} upper bound of 95% confidence interval{p_end}
{synopt:{cmd:r(threshold)}} optimal threshold (using interquartile range method)
{p_end}
{synopt:{cmd:r(pct_ge_thres)}} % of location pairs with correlations above the 
optimal threshold{p_end}
{synopt:{cmd:r(pct_ge_thres_nocl)}} % of inter-cluster location pairs with 
correlations above the optimal threshold{p_end}
{synopt:{cmd:r(T)}} number of time periods{p_end}
{synopt:{cmd:r(N)}} number of observations{p_end}
{synopt:{cmd:r(N_loc)}} number of locations{p_end}
{synopt:{cmd:r(N_clust)}} number of clusters{p_end}
{synopt:{cmd:r(N_outcomes)}} number of outcomes used to estimate correlations
between locations{p_end}
{synopt:{cmd:r(dof)}} degrees of freedom for estimating correlations between 
locations{p_end}
{synopt:{cmd:r(finite_sample_dof)}} finite sample adjustment for variance
calculation{p_end}
{p2colreset}{...}

{marker disclaimer}{...}
{title:Disclaimer}

{pstd}
This package is in beta/testing mode. Please use cautiously and feel free to 
report any errors to wjnkim@stanford.edu.

{pstd}
{cmd:tmo} requires the {cmd:{help gtools}} package. Please run {stata ssc install gtools} to install.

{marker references}{...}
{title:References}

{marker paper}{...}
{phang}DellaVigna, Stefano, Guido Imbens, Woojin Kim, and David Ritzwoller. 
(2025). "Using Multiple Outcomes to Adjust Standard Errors for Spatial 
Correlation."{p_end}
