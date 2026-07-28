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

* Load scpc FIRST, then tmo: scpc.ado runs -mata mata clear- when it loads,
* which would wipe tmo's mata functions (corr_resid etc., error r(3499)) if
* scpc were auto-loaded later, e.g. inside tmo's scpc_cmd() option.
qui do "$ROOT/scpc_tmo/scpc.ado"
qui do "$ADO/dev/tmo.ado"

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
*--- (0b) Distance cutoff shared by the Conley map (3) and the combined
*    TMO+Conley maps (6). Keep these equal: the panels of the combined figure
*    are read side by side, so a different radius in each is misleading.
global CONLEYCUT 150

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

* Plot the Conley Bartlett weight (decays with distance; 0 beyond the cutoff)
* rather than raw distance, so colour intensity = weight in the estimator.
* NB: $CONLEYCUT is used by BOTH this map and the combined TMO+Conley map of
* section 6, so that the two figures are directly comparable; it also matches
* the 150-mile bandwidth used in the applications of the paper.
gen conley_w = 1 - distance/$CONLEYCUT if distance<=$CONLEYCUT

#delimit ;
geoplot (area counties conley_w if `mainland', lwidth(none) lcolor(white)
            levels(6) missing(label("No weight")))
        (label counties NAME if `focal_cond', color(white) size(vsmall))
        (line states   if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/Conley.pdf", replace

*-------------------------------------------------------------------------------
*--- (4) One combined TMO+SCPC run feeds maps (4), (5) and (7)
*    Everything below comes from tmo's OWN internals, saved pair-level in the
*    dyad: corr (auxiliary-outcome correlations), Wfin (the SCPC kernel that
*    tmo computes internally when scpc_cmd() is given), and e(threshold).
*-------------------------------------------------------------------------------
preserve
use "$DAT/maps/cb_2018_us_county_20m.dta", clear
destring GEOID, replace
rename GEOID fips
keep fips _CX _CY
tempfile cent
save `cent'

use "$EXA/county_differences.dta", clear
merge 1:1 fips using `cent', keep(3) nogen
qui ds fips stfips PIN_persincpc_d EDU_college_d _CX _CY, not
local ylist `r(varlist)'

* scpc_cmd is incompatible with weights; savedyad stores the full dyad
qui tmo, cmd(reg PIN_persincpc_d EDU_college_d) x(EDU_college_d) ///
    ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    scpc_cmd(reg PIN_persincpc_d EDU_college_d, r) ///
    savedyad file("$TMP/mapsS")
scalar thresS = e(threshold)

use "$TMP/mapsS_dyad.dta", clear
keep if (id1==`focalfips' | id2==`focalfips') & id1!=id2
gen GEOID = cond(id1==`focalfips', id2, id1)
gen byte tmoselS = (abs(corr)>=thresS) & !missing(corr)
gen corrS = abs(corr)
gen aWfin = abs(Wfin)
keep GEOID tmoselS corrS aWfin
tempfile combS
save `combS'
restore

frame change counties
destring GEOID, replace
merge 1:1 GEOID using `combS', keep(1 3) nogen

*--- (4a) SCPC map: tmo's internal SCPC kernel weight |Wfin| for the focal row
xtile scpc_bin = aWfin, nq(6)

#delimit ;
geoplot (area counties scpc_bin if `mainland', lwidth(none) lcolor(white) levels(6)
        label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6") missing(label("No weight")))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr
graph export "$OUT/SCPC.pdf", replace

*-------------------------------------------------------------------------------
*--- (5) TMO map: only the pairs the estimator keeps (|corr| >= threshold)
*-------------------------------------------------------------------------------
gen tmo_kept = corrS if tmoselS==1

#delimit ;
geoplot (area counties tmo_kept if `mainland', lwidth(none) lcolor(white)
        levels(4) missing(label("No weight")))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo.pdf", replace

*--- (5b) Optional: rank counties by number of above-threshold partners
*    Useful for choosing the focal county of the maps. Set to 1 to run.
*    (2026-07: Tulsa=33 pairs; top counties: Winnebago WI=192, Walworth WI=191,
*     Jefferson WI=188, Wood OH=169 -- Great Lakes manufacturing belt)
local findfocal 0
if `findfocal'==1 {
    preserve
    use "$TMP/mapsS_dyad.dta", clear
    keep if id1!=id2
    gen byte above = (abs(corr)>=thresS) & !missing(corr)
    * stack so each pair counts for both of its counties
    expand 2, gen(dup)
    gen fips = cond(dup==0, id1, id2)
    collapse (sum) nabove=above, by(fips)
    gsort -nabove
    di as result _n "Top 10 counties by number of pairs above the TMO threshold:"
    li fips nabove in 1/10, noobs clean
    restore
}

*-------------------------------------------------------------------------------
*--- (6) Combined TMO + Conley from an ACTUAL combined run
*    tmo with distthreshold() computes pairwise distances internally (geodist)
*    and saves them in the dyad; the SE lets pairs enter either because
*    |corr| >= threshold (TMO) or because dist <= cutoff (Conley disk).
*-------------------------------------------------------------------------------
local conleycut $CONLEYCUT

preserve
use "$DAT/maps/cb_2018_us_county_20m.dta", clear
destring GEOID, replace
rename GEOID fips
keep fips _CX _CY
tempfile cent2
save `cent2'

use "$EXA/county_differences.dta", clear
merge 1:1 fips using `cent2', keep(3) nogen
qui ds fips stfips PIN_persincpc_d EDU_college_d _CX _CY, not
local ylist `r(varlist)'

qui tmo, cmd(reg PIN_persincpc_d EDU_college_d) x(EDU_college_d) ///
    ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    distthreshold(`conleycut') miles ///
    savedyad file("$TMP/mapsC")
scalar thresC = e(threshold)

use "$TMP/mapsC_dyad.dta", clear
keep if (id1==`focalfips' | id2==`focalfips') & id1!=id2
gen GEOID = cond(id1==`focalfips', id2, id1)
gen byte tmoselC = (abs(corr)>=thresC) & !missing(corr)
rename dist distC
keep GEOID tmoselC distC
tempfile combC
save `combC'
restore

merge 1:1 GEOID using `combC', keep(1 3) nogen

*--- (6a) uniform kernel (what distthreshold() does by default): categorical
gen cat = .
replace cat = 1 if distC<=`conleycut'
replace cat = 2 if tmoselC==1 & (distC>`conleycut' | missing(distC))
replace cat = 3 if tmoselC==1 & distC<=`conleycut'
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

*--- (6b) bartlett kernel (distkernel(bartlett)): non-selected pairs enter
*    with weight 1 - d/b, TMO pairs with full weight (shown in red)
* zero-weight pairs (beyond the cutoff) stay missing -> "No weight" grey
gen bartw = 1 - distC/`conleycut' if tmoselC!=1 & distC<`conleycut'

#delimit ;
geoplot (area counties bartw if `mainland' & tmoselC!=1, lwidth(none)
        lcolor(white) levels(6) missing(label("No weight")))
        (area counties if `mainland' & tmoselC==1, color(red) lwidth(none)
        label("TMO pair"))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo_conley_bartlett.pdf", replace

*-------------------------------------------------------------------------------
*--- (7) Combined TMO + SCPC map (same run as section 4)
*    TMO-selected pairs enter with their raw covariance (red); every other
*    pair enters through tmo's internal SCPC kernel weight Wfin (background).
*-------------------------------------------------------------------------------
xtile scpc_kbin = aWfin if tmoselS!=1, nq(6)

#delimit ;
geoplot (area counties scpc_kbin if `mainland' & tmoselS!=1, lwidth(none)
        lcolor(white) levels(6) label(1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6")
        missing(label("No weight")))
        (area counties if `mainland' & tmoselS==1, color(red) lwidth(none)
        label("TMO pair"))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/tmo_scpc.pdf", replace

*-------------------------------------------------------------------------------
*--- (8) Combined TMO + state clustering map
*    Pairs entering the variance: the focal county's own state (cluster) plus
*    the TMO-selected pairs from the same run as sections (4)-(5).
*-------------------------------------------------------------------------------
#delimit ;
geoplot (area counties if `mainland', color(gs14) lwidth(none)
        label("No weight"))
        (area counties if `mainland' & STATEFP=="`focalST'" & tmoselS!=1,
        color(ebblue) lwidth(none) label("Cluster (state)"))
        (area counties if `mainland' & tmoselS==1, color(red) lwidth(none)
        label("TMO pair"))
        (label counties NAME if `focal_cond', color(black) size(vsmall))
        (line states if `mainland', lwidth(vthin))
        , project legend(pos(5));
#delimit cr

graph export "$OUT/cluster_tmo.pdf", replace

* Clean up
frame change default
