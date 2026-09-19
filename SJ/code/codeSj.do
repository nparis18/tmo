/* codeSj.do                                                    tmo SJ examples
----|----1----|----2----|----3----|----4----|----5----|----6----|----7----|----8

  Generates ALL example logs (sjlog) and figures used in SJ/paper/spatial.tex
  Outputs:
    - Example logs  -> SJ/paper/examples  (files ending in .tex.log.tex,
      input by spatial.tex)
    - tmo figures   -> SJ/paper/figures/_hist.png and _qt.pdf
  Maps figures (state_clustering, Conley, SCPC, tmo) are made by maps.do
  NB: never write the sequence slash-star inside comments here -- Stata block
  comments nest, and an unbalanced open swallows the rest of the file.
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
global SJ   "$ROOT/SJ"
global DAT  "$ROOT/example"
global PPR  "$SJ/paper"
global FIG  "$SJ/paper/figures"
global TMP  "$SJ/temp"

* Which tmo version to use. Both files define a program called -tmo-, so the
* sjlog output always displays the command as tmo. NEVER load both in one
* session: their mata functions collide (e.g. corr_resid() signatures differ).
*   $ROOT/src      = release version (no areg/ivregress support)
*   $ROOT/src/dev  = development version (faster panel/IV; adds areg, ivregress)
global TMOSRC "$ROOT/src/dev"

* Load scpc FIRST, then tmo: scpc.ado runs -mata mata clear- when it loads,
* which wipes tmo's mata functions (error r(3499)) if loaded afterwards.
qui do "$ROOT/scpc_tmo/scpc.ado"
qui do "$TMOSRC/tmo.ado"

* Run from the paper folder so paths displayed inside the sjlogs are short
* and machine-independent (e.g. file("figures/"))
cd "$PPR"

*-------------------------------------------------------------------------------
*--- (1) County example (OLS): main illustrative example + figures
*-------------------------------------------------------------------------------
use "$DAT/county_differences.dta", clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

sjlog using "$PPR/examples/countyexample.tex", replace
tmo, cmd(reg PIN_persincpc_d EDU_college_d i.stfips, r) x(EDU_college_d) ///
    ylist(`ylist') i(fips) plothist plotq file("figures/")
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (2) Alternative estimation commands (regress, reghdfe, areg)
*-------------------------------------------------------------------------------
use "$DAT/county_differences.dta", clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

sjlog using "$PPR/examples/cmdexamples.tex", replace
timer clear
timer on 1
qui tmo, cmd(regress PIN_persincpc_d EDU_college_d i.stfips, r) ///
    x(EDU_college_d) ylist(`ylist') i(fips)
timer off 1

timer on 2
qui tmo, cmd(reghdfe PIN_persincpc_d EDU_college_d, vce(r) abs(stfips)) ///
    x(EDU_college_d) ylist(`ylist') i(fips)
timer off 2

timer on 3
qui tmo, cmd(areg PIN_persincpc_d EDU_college_d, r abs(stfips)) ///
    x(EDU_college_d) ylist(`ylist') i(fips)
timer off 3
timer list
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (3) Panel example
*-------------------------------------------------------------------------------
use "$DAT/county_panel.dta", clear
qui ds fips stfips EMN_farm EDU_publicenroll year, not
local ylist `r(varlist)'

sjlog using "$PPR/examples/panelexample.tex", replace
tmo, cmd(reg EMN_farm EDU_publicenroll i.year i.stfips, cluster(fips)) ///
    x(EDU_publicenroll) ylist(`ylist') i(fips) t(year)
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (4) IV example
*-------------------------------------------------------------------------------
use "$DAT/county_differences.dta", clear
qui ds fips stfips life_d VST_infmort_d AHRQ_emerdist_d AHRQ_obgyndist_d ///
    AHRQ_pediadist_d, not
local ylist `r(varlist)'

sjlog using "$PPR/examples/Ivexample.tex", replace
tmo, cmd(ivreg2 life_d (VST_infmort_d = AHRQ_emerdist_d AHRQ_obgyndist_d ///
    AHRQ_pediadist_d)) x(VST_infmort_d) ylist(`ylist') i(fips)
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (5) Combining tmo with other procedures
*-------------------------------------------------------------------------------
use "$DAT/county_differences.dta", clear
qui ds fips stfips PIN_persincpc_d EDU_college_d, not
local ylist `r(varlist)'

sjlog using "$PPR/examples/clusterExample.tex", replace
reg PIN_persincpc_d EDU_college_d, cluster(stfips)
tmo, cmd(reg PIN_persincpc_d EDU_college_d, cluster(stfips)) ///
    x(EDU_college_d) ylist(`ylist') i(fips)
sjlog close, replace

* Add county centroids for distance- and SCPC-based examples
preserve
use "$SJ/data/maps/cb_2018_us_county_20m.dta", clear
destring GEOID, replace
rename GEOID fips
keep fips _CX _CY
tempfile maps
save `maps'
restore

merge 1:1 fips using `maps', keep(3) nogen
qui ds fips stfips PIN_persincpc_d EDU_college_d _CX _CY, not
local ylist `r(varlist)'

* NB: the 150-mile bandwidth matches the maps of maps.do and the applications
* below, so that every distance-based figure and table in the paper uses the
* same cutoff
sjlog using "$PPR/examples/conleyExample.tex", replace
tmo, cmd(reg PIN_persincpc_d EDU_college_d, r) ///
    x(EDU_college_d) ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    distthreshold(150) miles
tmo, cmd(reg PIN_persincpc_d EDU_college_d, r) ///
    x(EDU_college_d) ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    distthreshold(150) miles distkernel(bartlett)
sjlog close, replace

// _CX: longitude ; _CY: latitude
rename (_CY _CX) (s_1 s_2)
// scpc reads s_* in physical variable order: force s_1 (lat) before s_2
order s_1 s_2

sjlog using "$PPR/examples/scpcExample.tex", replace
reg PIN_persincpc_d EDU_college_d, r
scpc, latlong
rename (s_1 s_2) (_CY _CX)

tmo, cmd(regress PIN_persincpc_d EDU_college_d, r) ///
    x(EDU_college_d) ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    scpc_cmd(reg PIN_persincpc_d EDU_college_d, r)
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (6) Second application: IV -- Bazzi et al. (2023), Southern white
*    migration and the vote for Trump. The dataset is rebuilt from the paper's
*    replication package (see ../data/replication/README_generated.txt).
*    Replicates Table 2 Panel A Column 4: shift-share IV, baseline controls,
*    state and truncation-dummy FEs, SEs clustered by 60x60-mile grid cell.
*
*    Every sjlog in sections (6) and (7) is printed in the paper next to the
*    table it produces, so each block must stay self-contained and readable,
*    and the numbers it displays must be exactly those in the table.
*-------------------------------------------------------------------------------

*--- (6a) Steps 1-2: baseline and auxiliary outcomes.
*    NOT printed in the paper: ordinary data preparation.
use "../data/replication/Republican_vote_data.dta", clear
keep if year==1940
local controls lnpopdens_hist pct_mfgempl pct_unemploy pct_laborforce ///
    pct_btot pct_popmexico pct_popgerman pct_popcanada pct_popireland ///
    pct_popitaly pct_farmacres pct_farmvalue pct_wilson_12 ///
    pct_cwenlistment pct_cwmortality
local iv (pct_Southerners_white = iv_mig_hat1_00_2)
local fe absorb(statefip dummy_trunmig_hat1_00_2)
local spec ivreghdfe Trump_share `iv' `controls', cluster(km_grid_cel_code) `fe'
quietly `spec'
scalar se0 = _se[pct_Southerners_white]
display "beta = " %5.3f _b[pct_Southerners_white] "   clustered SE = " %5.3f se0
generate byte insample = e(sample)
* candidates: all numeric variables except identifiers and geography, the
* outcome, the regressor, every ingredient of the shift-share instrument,
* missing-data flags and the controls; the loop then drops variables that are
* mostly missing, constant, or correlated above 0.8 with the model's variables
quietly ds icpsrfip_1 icpsrfip year fips county statefip South border North  ///
    D unincorp_1860 xcoord ycoord clon10 clat10 decade decade_1 area_sqmi*    ///
    km_grid_cel_code insample Trump_share votes_1940 votes_1948 votes_2000    ///
    votes_2016 tot_vote_pres candidatevotes totalvotes elec* Prep1900 FEpair* ///
    popdens_hist popdens_hist_00 Southerners_white* pct_Southerners_white*   ///
    Dpct_Southerners_white* D_pct_Southerners* iv_* dummy_trun* *_share_1900 ///
    Southerners_black_? Southerners_black_?? pct_Southerners_black_?          ///
    pct_Southerners_black_?? Northerners_white_?? pct_Northerners_white_??    ///
    Dpct_Northerners_white_?? pct_Northerners_white_hat* missing_* `controls', not
quietly ds `r(varlist)', has(type numeric)
local candidates `r(varlist)'
local ylist
foreach v of local candidates {
    quietly count if missing(`v') & insample
    if r(N) >= 0.5*e(N) continue
    quietly summarize `v' if insample
    if r(sd) == 0 continue
    local maxcorr 0
    foreach r in Trump_share pct_Southerners_white `controls' {
        quietly correlate `v' `r' if insample
        local maxcorr = max(`maxcorr', abs(r(rho)))
    }
    if `maxcorr' < 0.8 local ylist `ylist' `v'
}
display "auxiliary outcomes: " wordcount("`ylist'")

*--- (6b) Step 3: TMO augmenting the original grid clustering, with diagnostics
sjlog using "$PPR/examples/bazziExample.tex", replace
tmo, cmd(`spec') x(pct_Southerners_white) ylist(`ylist') i(fips) ///
    misslimit(0.5) plothist plothistnbins(100) plotse file("figures/bazzi")
sjlog close, replace
cap erase "$PPR/figures/bazzi_dyad.dta"

*--- (6c) Step 5: pair-level file behind Table bazzipredictors.  Only the tmo
*    call is printed; the merge and the tabulation are ordinary data work.
local opts x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5)
tempfile bzbase
quietly save `bzbase'

sjlog using "$PPR/examples/bazziPredictors.tex", replace
tmo, cmd(`spec') `opts' lat(clat10) lon(clon10) distthreshold(150) miles ///
    savedyad file("bazzi")
local threshold = e(threshold)
use bazzi_dyad.dta, clear
describe id1 id2 corr dist
sjlog close, replace

use `bzbase', clear
keep fips countyname statefip km_grid_cel_code totpop
rename (fips countyname statefip km_grid_cel_code totpop) (id1 name1 state1 grid1 pop1)
tempfile loc1 loc2
quietly save `loc1'
rename *1 *2
quietly save `loc2'
use bazzi_dyad.dta, clear
quietly drop if id1==id2
generate byte selected = abs(corr) >= `threshold' & !missing(corr)
quietly merge m:1 id1 using `loc1', keep(match) nogenerate
quietly merge m:1 id2 using `loc2', keep(match) nogenerate
generate byte within150 = dist <= 150
generate byte samegrid  = grid1 == grid2
generate byte samestate = state1 == state2
generate dpop = abs(pop1 - pop2)
quietly _pctile dpop, p(10)
generate byte closepop = dpop <= r(r1) if !missing(dpop)
generate byte any = within150 | samegrid | samestate | closepop
display "criterion" _col(14) "selected pairs" _col(32) "all pairs"
foreach v in within150 samegrid samestate closepop any {
    quietly summarize `v' if selected
    local sel = 100*r(mean)
    quietly summarize `v'
    display "`v'" _col(14) %4.1f `sel' "%" _col(32) %4.1f 100*r(mean) "%"
}
gsort -selected -dist
list name1 name2 corr dist in 1/3, noobs clean
use `bzbase', clear
cap erase "$PPR/bazzi_dyad.dta"

*--- (6d) Step 6: method comparison (Table bazzicompare)
sjlog using "$PPR/examples/bazziCompare.tex", replace
capture program drop tmoline
program tmoline
    args label option
    local se = cond("`option'"=="scpc", e(scpc_se), e(tmo_se))
    display "`label'" _col(16) "SE = " %5.3f `se' "   ratio = " %4.2f `se'/se0 _continue
    if "`option'"=="" display "   pairs = " %3.1f e(pct_ge_thres) "%"
    else display
end
local robust  ivreghdfe Trump_share `iv' `controls', robust `fe'
local bystate ivreghdfe Trump_share `iv' `controls', cluster(statefip) `fe'
quietly tmo, cmd(`robust') `opts'
tmoline "TMO"
quietly tmo, cmd(`spec') `opts'
tmoline "TMO + grid"
quietly tmo, cmd(`bystate') `opts'
tmoline "TMO + state"
quietly tmo, cmd(`spec') `opts' lat(clat10) lon(clon10) distthreshold(150) miles thresholdoff
tmoline "Conley"
quietly tmo, cmd(`spec') `opts' lat(clat10) lon(clon10) distthreshold(150) miles
tmoline "TMO + Conley"
quietly tmo, cmd(`spec') `opts' lat(clat10) lon(clon10) ///
    scpc_cmd(ivregress 2sls Trump_share `controls' i.statefip dummy_trunmig_hat1_00_2 `iv', robust)
tmoline "SCPC" scpc
tmoline "TMO + SCPC" nopairs
sjlog close, replace

*--- (6e) Step 7: sensitivity to the outcome collection (Table bazziylist)
* the trimmed collection drops the paper's secondary control sets and the vote
* counts; building it is data work and is not printed in the paper
quietly ds lnpopdens_hist_00 pct_*_00 share_breckinridge pct_bryan_96 oil_1900 ///
    oil_1940 AnyMines cottonmed potential_ag_prod d_coa d_riv d_lak elev_mean ///
    tri_ave tye_tfe890_500k_100_l6 votes_*
local secondary `r(varlist)'
local ytrim : list ylist - secondary

sjlog using "$PPR/examples/bazziYlist.tex", replace
* misuse: add measures of the regressor itself
local ymisuse `ylist' pct_Southerners_white1900 pct_Southerners_white_brdr ///
    Dpct_Southerners_white D_pct_Southerners_white_00_40
foreach collection in ylist ytrim ymisuse {
    quietly tmo, cmd(`spec') x(pct_Southerners_white) ylist(``collection'') ///
        i(fips) misslimit(0.5)
    display "`collection'" _col(12) "d = " wordcount("``collection''") ///
        "   df = " %4.1f e(dof) "   threshold = " %4.2f e(threshold) ///
        "   ratio = " %4.2f e(tmo_se)/se0
}
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (7) Step-by-step guide application: Bernini et al. (2023)
*    Rules for the auxiliary outcomes follow DellaVigna et al. (2025, Appendix
*    E.2); the 0.8 correlation cutoff is ours, as the rule gives no number.
*-------------------------------------------------------------------------------

*--- (7a) Steps 1-2: baseline and auxiliary outcomes.
*    NOT printed in the paper: this is ordinary data preparation, and the
*    guide shows only the code that uses the command itself.
use "../data/replication/dataset_wide_1.dta", clear
generate long fips = real(county)
local controls urbanB60 unemp60 family_less_3000 pop60 school_low ///
    cotton_suitability cotton_share_land1964 anti_black_county ///
    pro_black_county rep_share_1964
local spec reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.(`controls')#literacy_nc ibn.STATE, nocon robust cluster(judicial_divisions_id)
quietly `spec'
scalar se0 = _se[black_share60_lit_nc]
display "theta = " %5.3f _b[black_share60_lit_nc] "   clustered SE = " %5.3f se0
generate byte insample = e(sample)
* candidates: all numeric variables except identifiers, the outcome and
* regressor families (with their literacy interactions) and the controls
quietly ds county countycode FIPSTATE STATE geo judicial_divisions_id fips  ///
    literacy_nc SMD MIXED AL insample ch_ShareBl_* ShareBl_* diffshareblack* ///
    black_share* *_lit *_lit_* dist_lit* ln_rep_share_1964 `controls', not
quietly ds `r(varlist)', has(type numeric)
local candidates `r(varlist)'
local ylist
local ynocutoff
foreach v of local candidates {
    quietly count if missing(`v') & insample
    if r(N) >= 0.5*e(N) continue
    quietly summarize `v' if insample
    if r(sd) == 0 continue
    local ynocutoff `ynocutoff' `v'
    local maxcorr 0
    foreach r in ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 `controls' {
        quietly correlate `v' `r' if insample
        local maxcorr = max(`maxcorr', abs(r(rho)))
    }
    if `maxcorr' < 0.8 local ylist `ylist' `v'
}
display "auxiliary outcomes: " wordcount("`ylist'") ///
    "   (without the correlation cutoff: " wordcount("`ynocutoff'") ")"

*--- (7b) Step 3: TMO augmenting the original judicial-division clustering
sjlog using "$PPR/examples/berniniExample.tex", replace
tmo, cmd(`spec') x(black_share60_lit_nc) ylist(`ylist') i(fips) ///
    misslimit(0.5) plothist plothistnbins(100) plotse file("figures/bernini")
sjlog close, replace
cap erase "$PPR/figures/bernini_dyad.dta"

*--- (7c) Step 5: pair-level file behind Table bernpredictors.  Only the tmo
*    call is printed; merging county characteristics into the saved file and
*    tabulating the shares is ordinary data work.
preserve
use GEOID _CX _CY using "../data/maps/cb_2018_us_county_20m.dta", clear
generate long fips = real(GEOID)
tempfile centroids
quietly save `centroids'
restore
quietly merge 1:1 fips using `centroids', keep(master match) nogenerate
local opts x(black_share60_lit_nc) ylist(`ylist') i(fips) misslimit(0.5)
tempfile bernbase
quietly save `bernbase'

sjlog using "$PPR/examples/berniniPredictors.tex", replace
tmo, cmd(`spec') `opts' lat(_CY) lon(_CX) distthreshold(150) miles ///
    savedyad file("bernini")
local threshold = e(threshold)
use bernini_dyad.dta, clear
describe id1 id2 corr dist
sjlog close, replace

use `bernbase', clear
keep fips STATE judicial_divisions_id pop60 family_less_3000 urbanB60
rename (fips STATE judicial_divisions_id pop60 family_less_3000 urbanB60) ///
    (id1 state1 jdiv1 pop1 pov1 urb1)
tempfile loc1 loc2
quietly save `loc1'
rename *1 *2
quietly save `loc2'
use bernini_dyad.dta, clear
quietly drop if id1==id2
generate byte selected = abs(corr) >= `threshold' & !missing(corr)
quietly merge m:1 id1 using `loc1', keep(match) nogenerate
quietly merge m:1 id2 using `loc2', keep(match) nogenerate
generate byte within150 = dist <= 150
generate byte samestate = state1 == state2
generate byte samejdiv  = jdiv1 == jdiv2
foreach p in pop pov urb {
    generate d`p' = abs(`p'1 - `p'2)
    quietly _pctile d`p', p(10)
    generate byte close`p' = d`p' <= r(r1) if !missing(d`p')
}
generate byte any = within150 | samestate | samejdiv | closepop | closepov | closeurb
display "criterion" _col(14) "selected pairs" _col(32) "all pairs"
foreach v in within150 samestate samejdiv closepop closepov closeurb any {
    quietly summarize `v' if selected
    local sel = 100*r(mean)
    quietly summarize `v'
    display "`v'" _col(14) %4.1f `sel' "%" _col(32) %4.1f 100*r(mean) "%"
}
use `bernbase', clear
cap erase "$PPR/bernini_dyad.dta"

*--- (7d) Step 6: method comparison (Table berncompare)
sjlog using "$PPR/examples/berniniCompare.tex", replace
capture program drop tmoline
program tmoline
    args label option
    local se = cond("`option'"=="scpc", e(scpc_se), e(tmo_se))
    display "`label'" _col(16) "SE = " %5.3f `se' "   ratio = " %4.2f `se'/se0 _continue
    if "`option'"=="" display "   pairs = " %3.1f e(pct_ge_thres) "%"
    else display
end
quietly tmo, cmd(`spec') `opts' lat(_CY) lon(_CX) distthreshold(150) miles thresholdoff
tmoline "Conley"
quietly tmo, cmd(`spec') `opts'
tmoline "TMO"
quietly tmo, cmd(`spec') `opts' lat(_CY) lon(_CX) distthreshold(150) miles
tmoline "TMO + Conley"
sjlog close, replace

*--- (7e) Step 7: sensitivity to the outcome collection (Table bernylist)
sjlog using "$PPR/examples/berniniYlist.tex", replace
local ymisuse `ylist' black_share40 black_share50 diffshareblack60_50
foreach collection in ylist ynocutoff ymisuse {
    quietly tmo, cmd(`spec') x(black_share60_lit_nc) ylist(``collection'') ///
        i(fips) misslimit(0.5)
    display "`collection'" _col(12) "d = " wordcount("``collection''") ///
        "   df = " %4.1f e(dof) "   threshold = " %4.2f e(threshold) ///
        "   ratio = " %4.2f e(tmo_se)/se0
}
sjlog close, replace
