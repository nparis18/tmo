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
*global ROOT "/home/dcc213/code/tmo"
global ROOT "/Users/MacBook/Dropbox/Research/tmo_all/tmo"
global DAT "$ROOT/SJ/data"
global SRC "$ROOT/SJ/code"
global OUT "$ROOT/SJ/paper/figures"
global EXA "$ROOT/example"
global ADO "$ROOT/src"
global TMP "$ROOT/SJ/temp"

graph set window fontface "Arial Narrow"

local makeshape 0

* Focal county for ALL maps. Default: Tulsa, OK (fips 40143).
* Clearest TMO alternative: Winnebago, WI -> focalfips 55139, ST "55", CTY "139"
* (see the findfocal block in section 5b: Winnebago has 192 above-threshold
*  pairs vs Tulsa's 33). Change these three locals to re-centre every map.
local focalfips 40143
local focalST   "40"
local focalCTY  "143"
local focal_cond COUNTYFP=="`focalCTY'"&STATEFP=="`focalST'"

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
* Identify the focal county's state for highlighting
gen focal = 1 if STATEFP=="`focalST'"

frame change states
gen focal = 1 if STATEFP=="`focalST'"

#delimit ;
geoplot (area counties focal if `mainland', lwidth(none) lcolor(white) label("1") missing(label("No weight")))
        (label counties NAME if `focal_cond', color(white) size(vsmall))
        (line states   if `mainland', lwidth(vthin))
        , project ;
#delimit cr

graph export "$OUT/state_clustering.pdf", replace


*-------------------------------------------------------------------------------
*--- (3) Make Conley map
*-------------------------------------------------------------------------------
// spshape2dta creates centroid coordinates _CX and _CY
// We will use these to calculate distances
frame change counties
gen distance = .

// Set up the sp environment
spset, modify coordsys(latlong, miles)

* Find the _ID for the focal county
sum _ID if `focal_cond'
local focal_id = r(mean)

forvalues s=1/`=_N' {
    qui: spdistance `s' `focal_id'
    qui replace distance = r(distance) if _ID==`s'
}

* Plot the Conley Bartlett weight (decays with distance; 0 beyond 300 miles)
* rather than raw distance, so colour intensity = weight in the estimator
gen conley_w = 1 - distance/300 if distance<=300

#delimit ;
geoplot (area counties conley_w if `mainland', lwidth(none) lcolor(white)
            levels(6) missing(label("No weight")))
        (label counties NAME if `focal_cond', color(white) size(vsmall))
        (line states   if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/Conley.pdf", replace

*-------------------------------------------------------------------------------
*--- (4) Make SCPC map
*-------------------------------------------------------------------------------
preserve
use "$EXA/county_differences.dta", clear
rename fips GEOID

tempfile county_differences
save `county_differences'
restore

destring GEOID, replace
merge 1:1 GEOID using `county_differences', keep(3) nogen
* scpc latlong convention: s_1 = latitude, s_2 = longitude
* (spshape2dta centroids: _CX = longitude, _CY = latitude)
gen s_1 = _CY
gen s_2 = _CX

qui reg PIN_persincpc_d EDU_college_d, r
gen byte insample = e(sample)
scpc, latlong

* Wfin rows follow the estimation sample only; restrict accordingly
mata:
    W2 = Wfin[, 2..cols(Wfin)] // drop constant (col 1)
    P  = W2*W2'                // implicit SCPC pair-weight kernel

    sel = selectindex(st_data(., "insample"))
    id  = st_data(., "_ID", "insample")
    iFocal = selectindex(id:==`focal_id')[1]

    w = abs(P[iFocal, .])'
    st_addvar("double", "scpc_w")
    st_store(sel, "scpc_w", w)
end
xtile scpc_bin = scpc_w, nq(6)

#delimit ;
geoplot (area counties scpc_bin if `mainland', lwidth(none) lcolor(white) levels(6)
        label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6") missing(label("No weight")))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr
graph export "$OUT/SCPC.pdf", replace

*-------------------------------------------------------------------------------
*--- (5) Make TMO map
*-------------------------------------------------------------------------------
preserve
use "$EXA/county_differences.dta", clear

qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'
gen weight=1

qui do "$ADO/tmo.ado"
* NB: savedyad requires file(); dyad data saved as $TMP/maps_dyad.dta
qui tmo, cmd(reg PIN_persincpc_d EDU_college_d [fw = weight]) x(EDU_college_d) ylist(`ylist') i(fips) savedyad file("$TMP/maps")
scalar thres = e(threshold)
restore

*--- (5b) Optional: rank counties by number of above-threshold partners
*    Useful for choosing the focal county of the four maps. Set to 1 to run.
*    (2026-07: Tulsa=33 pairs; top counties: Winnebago WI=192, Walworth WI=191,
*     Jefferson WI=188, Wood OH=169 -- Great Lakes manufacturing belt)
local findfocal 0
if `findfocal'==1 {
    preserve
    use "$TMP/maps_dyad.dta", clear
    keep if id1!=id2
    gen byte above = (abs(corr)>=thres) & !missing(corr)
    * stack so each pair counts for both of its counties
    expand 2, gen(dup)
    gen fips = cond(dup==0, id1, id2)
    collapse (sum) nabove=above, by(fips)
    gsort -nabove
    di as result _n "Top 10 counties by number of pairs above the TMO threshold:"
    li fips nabove in 1/10, noobs clean
    restore
}

preserve
use "$TMP/maps_dyad.dta", clear
keep if (id2==`focalfips' | id1==`focalfips') & id1!=id2
gen GEOID = cond(id1==`focalfips', id2, id1)
keep GEOID corr
replace corr = abs(corr)

tempfile corr
save `corr'
restore

* Merge with counties frame
frame change counties
destring GEOID, replace
merge 1:1 GEOID using `corr', keep(1 3) nogen

* TMO keeps ONLY pairs with |corr| >= threshold: colour those, grey out the rest
gen byte tmosel = ((corr>=thres) & !missing(corr))
gen tmo_kept = corr if tmosel

#delimit ;
geoplot (area counties tmo_kept if `mainland', lwidth(none) lcolor(white)
        levels(4) missing(label("No weight")))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo.pdf", replace

*-------------------------------------------------------------------------------
*--- (6) Combined map: TMO + Conley (uniform kernel, `conleycut' miles)
*    Pairs entering the variance: Conley disk, TMO pairs beyond it, or both.
*    (The optimal threshold is unchanged when combining, so thres from (5)
*     applies; verified 2026-07: same .4329 with and without distthreshold.)
*-------------------------------------------------------------------------------
local conleycut 100

gen cat = .
replace cat = 1 if distance<=`conleycut'
replace cat = 2 if tmosel==1 & (distance>`conleycut' | missing(distance))
replace cat = 3 if tmosel==1 & distance<=`conleycut'
label define catl 1 "Conley (<=`conleycut' mi)" 2 "TMO pair" 3 "Both", replace
label values cat catl

#delimit ;
geoplot (area counties cat if `mainland', discrete lwidth(none) lcolor(white)
        missing(label("No weight")))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo_conley.pdf", replace

*-------------------------------------------------------------------------------
*--- (7) Combined map: TMO + SCPC
*    The combined estimator uses the raw covariance for TMO-selected pairs and
*    the SCPC kernel weight for all remaining pairs. Background: SCPC kernel
*    bins from (4); red: TMO-selected pairs from (5).
*-------------------------------------------------------------------------------
#delimit ;
geoplot (area counties scpc_bin if `mainland' & tmosel!=1, lwidth(none)
        lcolor(white) levels(6) label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6")
        missing(label("No weight")))
        (area counties if `mainland' & tmosel==1, color(red) lwidth(none)
        label("TMO pair"))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo_scpc.pdf", replace

* Clean up
frame change default
