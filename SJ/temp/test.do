clear all
cd "/Users/MacBook/Dropbox/Research/tmo"

*use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_differences.dta", clear
*save county_differences, replace

use county_differences, clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'
gen weight=1

qui reg PIN_persincpc_d EDU_college_d [fw = weight], r cluster(stfips)

*qui do "/Users/MacBook/Dropbox/Research/tmo/tmo.ado"
*tmo, cmd(reg PIN_persincpc_d EDU_college_d [fw = weight], r) x(EDU_college_d) ylist(`ylist') i(fips)

qui do "/Users/MacBook/Dropbox/Research/tmo/tmo_NP.ado"
timer clear
local avr1
forval i =1/100{
timer on 1
qui tmo_NP, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ylist(`ylist') i(fips)
timer off 1
}
quietly timer list 1      // deja en r(): r(t1) (total, seg) y r(nt1) (# de veces)
local prom = r(t1) / r(nt1)
display as text "Total (s): " %9.4f r(t1) ///
        "   #iteraciones: " r(nt1) ///
        "   Promedio (s): " %9.6f `prom'

        stop
//use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_panel.dta", clear
//save county_differences_panel, replace
/*
use county_differences_panel, clear
*keep if inrange(stfips, 1,10)

qui ds fips stfips EMN_farm EDU_publicenroll year, not
local ylist `r(varlist)'

do "/Users/MacBook/Dropbox/Research/tmo/tmo_NP.ado"
forval i=1/1 {
    timer clear
    timer on 1
    //log using "/Users/MacBook/Dropbox/Research/tmo_panel", t replace
    tmo_NP, cmd(areg EMN_farm EDU_publicenroll, absorb(stfips) cluster(fips)) ///
    x(EDU_publicenroll) ylist(`ylist') i(fips) t(year)
    //log close
    timer off 1
    timer list
}
*/
/*s
timer clear
timer on 1
tmo_NP, cmd(reghdfe EMN_farm EDU_publicenroll i.year, absorb(stfips) cluster(fips)) ///
x(EDU_publicenroll) ylist(`ylist') i(fips) t(year)
timer off 1

qui do "/Users/MacBook/Dropbox/Research/tmo/tmo.ado"
timer on 2
tmo, cmd(reghdfe EMN_farm EDU_publicenroll i.year, absorb(stfips) cluster(fips)) ///
x(EDU_publicenroll) ylist(`ylist') i(fips) t(year)
timer off 2

timer list
*/

use county_differences, clear
qui ds fips stfips life_d VST_infmort_d AHRQ_emerdist_d AHRQ_obgyndist_d AHRQ_pediadist_d, not
local ylist `r(varlist)'

tmo, cmd(ivreg2 life_d (VST_infmort_d = AHRQ_emerdist_d AHRQ_obgyndist_d AHRQ_pediadist_d)) x(VST_infmort_d) ylist(`ylist') i(fips)
