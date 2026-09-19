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
*    Replicates Table 2 Panel A Column 4.
*
*    As in section (7), the printed blocks are written out in full; the
*    selection code that produces the outcome list is not printed.
*-------------------------------------------------------------------------------

*--- (6a) data preparation (not printed)
use "../data/replication/Republican_vote_data.dta", clear
keep if year==1940
local controls lnpopdens_hist pct_mfgempl pct_unemploy pct_laborforce ///
    pct_btot pct_popmexico pct_popgerman pct_popcanada pct_popireland ///
    pct_popitaly pct_farmacres pct_farmvalue pct_wilson_12 ///
    pct_cwenlistment pct_cwmortality
quietly ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `controls', cluster(km_grid_cel_code) absorb(statefip dummy_trunmig_hat1_00_2)
generate byte insample = e(sample)
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
tempfile bzbase
quietly save `bzbase'

*--- (6b) Steps 1-3: the specification, the outcome list, and the first run
sjlog using "$PPR/examples/bazziExample.tex", replace
* the historical controls of Table 2, Panel A
global controls lnpopdens_hist pct_mfgempl pct_unemploy pct_laborforce ///
    pct_btot pct_popmexico pct_popgerman pct_popcanada pct_popireland ///
    pct_popitaly pct_farmacres pct_farmvalue pct_wilson_12 ///
    pct_cwenlistment pct_cwmortality

* the 85 auxiliary outcomes selected in Step 2
global outcomes Prep_28 Prep_32_64 Prep_72_00 Pdem_28 wallace_1968 swing_48_00 ///
    delta_Prep_1940_00 Pprog_12 Pprog_24 Pprog_48 Pdixie_68 pop natives ///
    natives_white natives_black Southerners_black pct_Southerners_black ///
    totpop wtot wnat wfor btot popcanada popmexico popireland popgerman ///
    popitaly farmempl farmnum farmvalue farmacres farmoutput farmten ///
    farmlarge farmblack mfgempl mfgwage mfgoutput unemploy laborforce ///
    votes_dem_92 votes_total_92 pct_cleveland_92 votes_dem_96 ///
    votes_total_96 pct_bryan_96 votes_dem_12 votes_total_12 urban ///
    elev_mean d_coa d_riv d_lak tri_ave potential_ag_prod ///
    tye_tfe890_500k_100_l6 pct_mfgempl_00 pct_popcanada_00 pct_popitaly_00 ///
    pct_farmacres_00 pct_farmvalue_00 alfalfa_suit wheat_suit pulses_suit ///
    cotton_suit potato_suit sweetpotato_suit oats_suit maize_suit ///
    tobacco_suit cottonmed votes_breckinridge share_breckinridge ///
    popmale1344 nbenlisted nbdead china_shock oil_1940 oil_1900 AnyMines ///
    pct_Southerners_black1900 popsqmi1860 Northerners_white ///
    pct_Northerners_white Dpct_Northerners_white

* the published specification: shift-share IV, state and truncation-dummy
* fixed effects, standard errors clustered by 60x60-mile grid cell
quietly ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    $controls, cluster(km_grid_cel_code) ///
    absorb(statefip dummy_trunmig_hat1_00_2)
display "coefficient " %5.3f _b[pct_Southerners_white] ///
    "    clustered SE " %5.3f _se[pct_Southerners_white]
scalar clustered = _se[pct_Southerners_white]

* the same specification with the TMO adjustment; cmd() takes the IV command
* unchanged and x() names the endogenous regressor
tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5) ///
    plothist plothistnbins(100) plotse file("figures/bazzi")
sjlog close, replace
cap erase "$PPR/figures/bazzi_dyad.dta"
local a $outcomes
local b `ylist'
local d1 : list a - b
local d2 : list b - a
assert "`d1'`d2'" == ""

*--- (6c) Step 5: the pair-level file behind Table bazzipredictors
use `bzbase', clear
sjlog using "$PPR/examples/bazziPredictors.tex", replace
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) distthreshold(150) miles savedyad file("bazzi")
display "selection threshold " %4.2f e(threshold)
use bazzi_dyad.dta, clear
describe id1 id2 corr dist
sjlog close, replace
local threshold = 0.51

* the merge and the tabulation are ordinary data work (not printed)
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

*--- (6d) Step 6: the corrections compared in Table bazzicompare
sjlog using "$PPR/examples/bazziCompare.tex", replace
* TMO on its own: the command passed to cmd() carries no clustering
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, robust ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5)
display "TMO            SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"

* TMO augmenting the original grid clustering
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5)
display "TMO + grid     SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"

* TMO augmenting state clusters: only the variance option changes
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(statefip) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5)
display "TMO + state    SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"

* Conley alone, and Conley combined with TMO
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) distthreshold(150) miles thresholdoff
display "Conley         SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) distthreshold(150) miles
display "TMO + Conley   SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"

* SCPC: scpc_cmd() takes the regression written for ivregress, which scpc
* supports; e(scpc_se) is SCPC on its own and e(tmo_se) the combination
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) ///
    scpc_cmd(ivregress 2sls Trump_share $controls i.statefip ///
    dummy_trunmig_hat1_00_2 (pct_Southerners_white = iv_mig_hat1_00_2), robust)
display "SCPC           SE " %5.3f e(scpc_se) "   ratio " %4.2f e(scpc_se)/clustered
display "TMO + SCPC     SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered
sjlog close, replace

*--- (6e) Step 7: the collections compared in Table bazziylist
* the trimmed collection drops the 1900 controls, the sorting confounds and
* the vote counts from $outcomes; building it is data work (not printed)
quietly ds lnpopdens_hist_00 pct_*_00 share_breckinridge pct_bryan_96 oil_1900 ///
    oil_1940 AnyMines cottonmed potential_ag_prod d_coa d_riv d_lak elev_mean ///
    tri_ave tye_tfe890_500k_100_l6 votes_*
local secondary `r(varlist)'
local ytrim : list ylist - secondary
global trimmed `ytrim'

sjlog using "$PPR/examples/bazziYlist.tex", replace
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($outcomes) i(fips) misslimit(0.5)
display "baseline           outcomes " e(N_outcomes)-1 "   df " %4.1f e(dof) ///
    "   threshold " %4.2f e(threshold) "   ratio " %4.2f e(tmo_se)/clustered

* $trimmed drops the 1900 controls, the sorting confounds and the vote counts
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist($trimmed) i(fips) misslimit(0.5)
display "trimmed            outcomes " e(N_outcomes)-1 "   df " %4.1f e(dof) ///
    "   threshold " %4.2f e(threshold) "   ratio " %4.2f e(tmo_se)/clustered

* misuse: four measures of the regressor itself added to the collection
quietly tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
        $controls, cluster(km_grid_cel_code) ///
        absorb(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ///
    ylist($outcomes pct_Southerners_white1900 pct_Southerners_white_brdr ///
    Dpct_Southerners_white D_pct_Southerners_white_00_40) ///
    i(fips) misslimit(0.5)
display "misuse             outcomes " e(N_outcomes)-1 "   df " %4.1f e(dof) ///
    "   threshold " %4.2f e(threshold) "   ratio " %4.2f e(tmo_se)/clustered
sjlog close, replace

*-------------------------------------------------------------------------------
*--- (7) Step-by-step guide application: Bernini et al. (2023)
*    Rules for the auxiliary outcomes follow DellaVigna et al. (2025, Appendix
*    E.2); the 0.8 correlation cutoff is ours, as the rule gives no number.
*
*    The blocks printed in the paper are written out in full, with the actual
*    variable names, so that a reader can copy them. The selection code that
*    produced the outcome list is below but is NOT printed: it is data work,
*    specific to this replication package. The -assert- lines guarantee that
*    the lists written out below are the ones those rules produce.
*-------------------------------------------------------------------------------

*--- (7a) data preparation (not printed)
use "../data/replication/dataset_wide_1.dta", clear
generate long fips = real(county)
local controls urbanB60 unemp60 family_less_3000 pop60 school_low ///
    cotton_suitability cotton_share_land1964 anti_black_county ///
    pro_black_county rep_share_1964
quietly reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.(`controls')#literacy_nc ibn.STATE, nocon robust cluster(judicial_divisions_id)
generate byte insample = e(sample)
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

* county centroids, for the distance-based runs
preserve
use GEOID _CX _CY using "../data/maps/cb_2018_us_county_20m.dta", clear
generate long fips = real(GEOID)
tempfile centroids
quietly save `centroids'
restore
quietly merge 1:1 fips using `centroids', keep(master match) nogenerate
tempfile bernbase
quietly save `bernbase'

*--- (7b) Steps 1-3: the specification, the outcome list, and the first run
sjlog using "$PPR/examples/berniniExample.tex", replace
* controls of Table 2, Column 4, each interacted with federal coverage
global controls urbanB60 unemp60 family_less_3000 pop60 school_low ///
    cotton_suitability cotton_share_land1964 anti_black_county ///
    pro_black_county rep_share_1964

* the 60 auxiliary outcomes selected in Step 2
global outcomes population60 all_officials_1964 ln_tnt_pres_tot_1940 ///
    ln_tnt_pres_tot_1960 county_expenditure57_re_pc county_current57_re_pc ///
    county_cap57_re_pc county_expenditure82_re_pc county_current82_re_pc ///
    county_cap82_re_pc county_cap_pre county_cap_post county_current_pre ///
    county_current_post county_tot_pre county_tot_post ///
    ch_AllCountyOfficials ch_AllMunicipality ch_AllEducation ///
    naacp_branch_1964_pcb naacp_branch_1942_pcb ch_naacp_64_42_pcb ///
    kkk_klavern_64_66_pcw kkk_klavern_15_40_pcw ch_kkk_66_40_pcw ///
    black_lynching_1930_1940 black_lynching_1950_1964 ///
    ch_black_lynching60_40 ln_tnt_gov_1940 ln_tnt_gov_1960 ///
    lndiff_tnt_gov_60_40 lndiff_pop_50_60 lndiff_pop_60_80 ///
    urbanchange60_50 urbanchange unemp50 unempchange ruralchange ///
    ruralchange60_50 school_lowchange ch_cotton_sh_land_64_45 ///
    lndiff_tnt_pres_60_40 ln_rep_share_1960 ln_rep_share_1952 ///
    ln_rep_share_1940 lndiff_rep_share_64_40 lndiff_rep_share_60_40 ///
    ch_black_share_60_50 lndiff_lower_60_50 lndiff_upper_60_50 ln_lower_60 ///
    ln_lower_50 ln_upper_60 ln_upper_50 unempchange60_50 ///
    familypovchange60_50 school_lowchange60_50 ln_winner_gov_60 ///
    ln_winner_gov_40 lndiff_winner_gov

* the published specification, clustered by judicial division
quietly reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.($controls)#literacy_nc ibn.STATE, nocon robust ///
    cluster(judicial_divisions_id)
display "coefficient " %5.3f _b[black_share60_lit_nc] ///
    "    clustered SE " %5.3f _se[black_share60_lit_nc]
scalar clustered = _se[black_share60_lit_nc]

* the same specification, with the TMO adjustment and the two diagnostic plots
tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes) i(fips) misslimit(0.5) ///
    plothist plothistnbins(100) plotse file("figures/bernini")
sjlog close, replace
cap erase "$PPR/figures/bernini_dyad.dta"
local a $outcomes
local b `ylist'
local d1 : list a - b
local d2 : list b - a
assert "`d1'`d2'" == ""

*--- (7c) Step 5: the pair-level file behind Table bernpredictors
use `bernbase', clear
sjlog using "$PPR/examples/berniniPredictors.tex", replace
* latitude(), longitude() and distthreshold() add the distance between the two
* locations to the file that savedyad writes, one observation per pair
quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(_CY) lon(_CX) distthreshold(150) miles savedyad file("bernini")
display "selection threshold " %4.2f e(threshold)
use bernini_dyad.dta, clear
describe id1 id2 corr dist
sjlog close, replace
local threshold = 0.52

* merging county characteristics into that file and tabulating the shares is
* ordinary data work and is not printed in the paper
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

*--- (7d) Step 6: the corrections compared in Table berncompare
sjlog using "$PPR/examples/berniniCompare.tex", replace
* Conley alone: thresholdoff switches the TMO selection off, leaving the
* distance band
quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(_CY) lon(_CX) distthreshold(150) miles thresholdoff
display "Conley         SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"

* TMO augmenting the original clustering
quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes) i(fips) misslimit(0.5)
display "TMO            SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"

* both: a pair enters if it is within the band or if its outcomes co-move
quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes) i(fips) misslimit(0.5) ///
    lat(_CY) lon(_CX) distthreshold(150) miles
display "TMO + Conley   SE " %5.3f e(tmo_se) "   ratio " %4.2f e(tmo_se)/clustered ///
    "   pairs " %4.1f e(pct_ge_thres) "%"
sjlog close, replace

*--- (7e) Step 7: the collections compared in Table bernylist
sjlog using "$PPR/examples/berniniYlist.tex", replace
* the seven outcomes that the 0.8 correlation cutoff removed
global dropped all_officials_1980 pop50 urban50 family_less50_2000 familypovchange ///
    school_low50 cotton_share_land1945

quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes) i(fips) misslimit(0.5)
display "baseline           outcomes " e(N_outcomes)-1 "   df " %4.1f e(dof) ///
    "   threshold " %4.2f e(threshold) "   ratio " %4.2f e(tmo_se)/clustered

quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist($outcomes $dropped) i(fips) misslimit(0.5)
display "without cutoff     outcomes " e(N_outcomes)-1 "   df " %4.1f e(dof) ///
    "   threshold " %4.2f e(threshold) "   ratio " %4.2f e(tmo_se)/clustered

* misuse: three measures of the regressor itself added to the collection
quietly tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.($controls)#literacy_nc ibn.STATE, nocon robust ///
        cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ///
    ylist($outcomes black_share40 black_share50 diffshareblack60_50) ///
    i(fips) misslimit(0.5)
display "misuse             outcomes " e(N_outcomes)-1 "   df " %4.1f e(dof) ///
    "   threshold " %4.2f e(threshold) "   ratio " %4.2f e(tmo_se)/clustered
sjlog close, replace
local a $outcomes $dropped
local b `ynocutoff'
local d1 : list a - b
local d2 : list b - a
assert "`d1'`d2'" == ""
