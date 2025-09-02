
global ROOT "/Users/MacBook/Dropbox/Research/tmo_all/tmo"
global DAT "$ROOT/example"
global SRC "$ROOT/src"

use "${DAT}/county_differences.dta", clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'
gen weight=1
 log using "/Users/MacBook/Dropbox/Research/tmo_all/tmo/SJ/temp/log/tempLog", t replace

qui do "${SRC}/tmo.ado"
tmo, cmd(reg PIN_persincpc_d EDU_college_d, cluster(stfips)) x(EDU_college_d) ylist(`ylist') i(fips) savedyad
log close

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




use "${DAT}/county_differences.dta", clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

local avr1
qui do "${SRC}/tmo.ado"
timer clear

forval i =1/100{
timer on 1
tmo, cmd(reg PIN_persincpc_d EDU_college_d, cluster(stfips)) x(EDU_college_d) ylist(`ylist') i(fips)
timer off 1
}
quietly timer list 1      // deja en r(): r(t1) (total, seg) y r(nt1) (# de veces)
local prom = r(t1) / r(nt1)
display as text "Total (s): " %9.4f r(t1) ///
        "   #iteraciones: " r(nt1) ///
        "   Promedio (s): " %9.6f `prom'

// 7.600730
use "${DAT}/county_panel.dta", clear
*keep if inrange(stfips, 1,10)

qui ds fips stfips EMN_farm EDU_publicenroll year, not
local ylist `r(varlist)'
qui do "${SRC}/tmo_new.ado"

tmo_NP, cmd(reghdfe EMN_farm EDU_publicenroll, absorb(stfips) cluster(stfips)) ///
x(EDU_publicenroll) ylist(`ylist') i(fips) t(year) plotq plothist plotse file(./example) threshold(0.5)

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



/*
1) I need to review if all the options works
    - Personalised threshold
    - IV (and panel IV)