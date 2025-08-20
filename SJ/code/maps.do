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

graph set window fontface "Arial Narrow"

local makeshape 0

*-------------------------------------------------------------------------------
*--- (1) Make shape files
*-------------------------------------------------------------------------------
if `makeshape'==1 {
    cd "$DAT/maps"
    spshape2dta "$DAT/maps/cb_2018_us_county_20m/cb_2018_us_county_20m", replace
    spshape2dta "$DAT/maps/cb_2018_us_state_20m/cb_2018_us_state_20m", replace
    cd "$SRC"
}
* Load county and state shape data into frames
geoframe create counties "$DAT/maps/cb_2018_us_county_20m", shp("$DAT/maps/cb_2018_us_county_20m_shp")
geoframe create states   "$DAT/maps/cb_2018_us_state_20m", shp("$DAT/maps/cb_2018_us_state_20m_shp")


*-------------------------------------------------------------------------------
*--- (2) Make clustering map
*-------------------------------------------------------------------------------

local mainland STATEFP!="02"&STATEFP!="15"&STATEFP!="72"

frame change counties
* Identify Tulsa county, Oklahoma, for highlighting
gen Tulsa = 1 if STATEFP=="40"
local Tulsa_cond COUNTYFP=="143"&STATEFP=="40"

frame change states
gen Tulsa = 1 if STATEFP=="40"

#delimit ;
geoplot (area counties Tulsa if `mainland', lwidth(none) lcolor(white) label("1") missing(label("No weight")))
        (label counties NAME if `Tulsa_cond', color(white) size(vsmall))
        (line states   if `mainland', lwidth(vthin))
        , project ;
#delimit cr

graph export "$OUT/state_clustering.pdf", replace


*-------------------------------------------------------------------------------
*--- (3) Make Conley map
*-------------------------------------------------------------------------------
* Tulsa coordinates
local ref_long = -95.940139
local ref_lat  = 36.119398

// spshape2dta creates centroid coordinates _CX and _CY
// We will use these to calculate distances
frame change counties
gen distance = .

// Set up the sp anvironment
spset, modify coordsys(latlong, miles)

* Find the _ID for Tulsa County, OK
sum _ID if `Tulsa_cond'
local tulsa_id = r(mean)

forvalues s=1/`=_N' {
    qui: spdistance `s' `tulsa_id'
    qui replace distance = r(distance) if _ID==`s'
}
replace distance = . if distance>300

#delimit ;
geoplot (area counties distance if `mainland', lwidth(none) lcolor(white) 
            label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6") missing(label("No weight")))
        (label counties NAME if `Tulsa_cond', color(white) size(vsmall))
        (line states   if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/Conley.pdf", replace

*-------------------------------------------------------------------------------
*--- (4) Make SCPC map
*-------------------------------------------------------------------------------
preserve
use county_differences, clear
rename fips GEOID

tempfile county_differences
save `county_differences'
restore

destring GEOID, replace
merge 1:1 GEOID using `county_differences', keep(3) nogen
gen s_1 = _CX
gen s_2 = _CY

qui reg PIN_persincpc_d EDU_college_d, r
scpc, latlong

mata:
    W2 = Wfin[, 2..cols(Wfin)] // drop constants
    P  = W2*W2'                // kernel SCPC

    id = st_data(., "_ID")
    iTulsa = selectindex(id:==`tulsa_id')[1]

    w = abs(P[iTulsa, .])'
    st_addvar("double", "scpc_w")
    st_store(., "scpc_w", w)
end
xtile scpc_bin = scpc_w, nq(6)

#delimit ;
geoplot (area counties scpc_bin if `mainland', lwidth(none) lcolor(white) levels(6)
        label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6") missing(label("No weight")))
        (label counties NAME if `Tulsa_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr
graph export "$OUT/SCPC.pdf", replace

*-------------------------------------------------------------------------------
*--- (5) Make TMO map
*-------------------------------------------------------------------------------
preserve
use "https://raw.githubusercontent.com/wjnkim/tmo/master/example/county_differences.dta", clear

qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'
gen weight=1

qui do "/Users/MacBook/Dropbox/Research/tmo/tmo.ado"
qui tmo, cmd(reg PIN_persincpc_d EDU_college_d [fw = weight]) x(EDU_college_d) ylist(`ylist') i(fips) savedyad
scalar thres = e(threshold)
restore

preserve
use "/Users/MacBook/Dropbox/Research/tmo/_dyad.dta", clear
keep if id2==40143 | id1==40143
gen idAux = id1 if id1>40143
replace id1=id2 if id1>40143
replace id2=idAux if idAux!=.
keep id2 corr
rename id2 GEOID
replace corr =abs(corr)
gen threshold= corr>thres

tempfile corr
save `corr'
restore

local mainland STATEFP!="02"&STATEFP!="15"&STATEFP!="72"
local Tulsa_cond COUNTYFP=="143"&STATEFP=="40"

* Merge with counties frame
frame change counties
destring GEOID, replace
merge 1:1 GEOID using `corr'

gen mark= ((abs(corr)>=thres) & !missing(corr))
xtile tmo_bin = corr, nq(6)

* Plot the TMO correlations
#delimit ;
geoplot (area counties tmo_bin if `mainland', lwidth(none) lcolor(white) levels(6)
        label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6") missing(label("No weight")))
        (label counties NAME if `Tulsa_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo.pdf", replace

* Clean up
drop COUNTYFP_str corr corr_category
frame change default
frame drop tmo_data
