/* maps.do                            DCC                  yyyy-mm-dd:2025-05-18
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

  Sets up maps to simply display the difference between methods
*/

vers 16
clear all
set more off
cap log close

*-------------------------------------------------------------------------------
*--- (0) Set up the environment
*-------------------------------------------------------------------------------
*global ROOT "/home/dcc213/code/tmo/SJ"
global ROOT "/Users/MacBook/Dropbox/Research/tmo/SJ"
global DAT "$ROOT/data"
global SRC "$ROOT/code"
global OUT "$ROOT/results/figures"

clear all
do "/Users/MacBook/Dropbox/Research/tmo/tmo_NP.ado"
*use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_differences.dta", clear
*save county_differences, replace
use county_differences, clear

qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

sjlog using "$SJ/cmdexamples.tex", replace
noisily timer on 1
qui tmo_NP, cmd(regress PIN_persincpc_d EDU_college_d i.stfips, r) /// 
x(EDU_college_d) ylist(`ylist') i(fips)
timer off 1

timer on 2
qui tmo_NP, cmd(reghdfe PIN_persincpc_d EDU_college_d, vce(r) abs(stfips)) /// 
x(EDU_college_d) ylist(`ylist') i(fips)
timer off 2

timer on 3
qui tmo_NP, cmd(areg PIN_persincpc_d EDU_college_d, r abs(stfips)) /// 
x(EDU_college_d) ylist(`ylist') i(fips)
timer off 3
timer list
sjlog close, replace

**************************************** Panel example*************************************************+

*sjlog using "$SJ/panelexample.tex", replace
use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_panel.dta", clear
qui do "/Users/MacBook/Dropbox/Research/tmo/tmo_NP.ado"

qui ds fips stfips EMN_farm EDU_publicenroll year, not
local ylist `r(varlist)'

tmo_NP, cmd(reg EMN_farm EDU_publicenroll i.year i.stfips , ///
cluster(fips)) x(EDU_publicenroll) ylist(`ylist') i(fips) t(year)
*sjlog close, replace

*************************************** IV example *****************************************************************

sjlog using "$SJ/Ivexample.tex", replace
use county_differences, clear
qui ds fips stfips life_d VST_infmort_d AHRQ_emerdist_d AHRQ_obgyndist_d AHRQ_pediadist_d, not
local ylist `r(varlist)'

tmo, cmd(ivreg2 life_d (VST_infmort_d = AHRQ_emerdist_d AHRQ_obgyndist_d AHRQ_pediadist_d)) ///
x(VST_infmort_d) ylist(`ylist') i(fips)
sjlog close, replace

********************************* Alternative estimators ***********************************************************

sjlog using "$SJ/alternativeExample.tex", replace
use county_differences, clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
global ylist `r(varlist)'

reg PIN_persincpc_d EDU_college_d, cluster(stfips)
tmo_NP, cmd(reg PIN_persincpc_d EDU_college_d, cluster(stfips))

preserve
use "$DAT/maps/cb_2018_us_county_20m.dta", clear
destring GEOID, replace

tempfile maps
save `maps'
restore

rename fips GEOID
merge 1:1 GEOID using `maps', keep(3) nogen
// longitude x
// latitude y
rename(_CY _CX)(s_1 s_2)
reg PIN_persincpc_d EDU_college_d, r
scpc, latlong
rename(s_1 s_2)(_CY _CX)

tmo, cmd(regress PIN_persincpc_d EDU_college_d, r) /// 
x(EDU_college_d) ylist(${ylist}) i(GEOID) lat(_CY) lon(_CX) scpc_cmd(reg PIN_persincpc_d EDU_college_d,r)

sjlog close, replace

/*
sjlog using "$SJ/countyexample.tex", replace
use ../example/county_differences
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips)
sjlog close, replace
exit


use ../example/county_differences
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

timer on 1
tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips)
timer off 1 
exit

timer on 2
tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips)
timer off 2

timer list



*/