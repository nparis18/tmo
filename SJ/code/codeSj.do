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

sjlog using "$PPR/examples/conleyExample.tex", replace
tmo, cmd(reg PIN_persincpc_d EDU_college_d, r) ///
    x(EDU_college_d) ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    distthreshold(100) miles
tmo, cmd(reg PIN_persincpc_d EDU_college_d, r) ///
    x(EDU_college_d) ylist(`ylist') i(fips) lat(_CY) lon(_CX) ///
    distthreshold(100) miles distkernel(bartlett)
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
*    migration and the vote for Trump. Fully self-contained: the dataset is
*    built from the paper's replication package (see README_generated.txt next
*    to the .dta) and all auxiliary outcomes come from that package.
*    Replicates Table 2 Panel A Column 4 (the main specification, per the
*    companion paper): shift-share IV, baseline controls, state and
*    truncation-dummy FEs, SEs clustered by 60x60-mile grid cell.
*    Published values: coef 1.03, cluster SE 0.17, n = 1,886.
*-------------------------------------------------------------------------------
use "$SJ/data/replication/Republican_vote_data.dta", clear
keep if year==1940

local ctrl lnpopdens_hist pct_mfgempl pct_unemploy pct_laborforce pct_btot ///
    pct_popmexico pct_popgerman pct_popcanada pct_popireland pct_popitaly ///
    pct_farmacres pct_farmvalue pct_wilson_12 pct_cwenlistment pct_cwmortality

* baseline: same estimates as BFFPT's Table 2 Panel A Column 4
qui ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) `ctrl', ///
    cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)
scalar se_bz_base = _se[pct_Southerners_white]
di as text "BFFPT Table 2 Panel A Col 4:  beta = " %6.4f _b[pct_Southerners_white] ///
    "   cluster SE = " %6.4f _se[pct_Southerners_white] "   N = " e(N)
gen byte insamp = e(sample)

* Step 2: auxiliary outcomes, applying the same rules as the guide: all
* feasible package variables, excluding (i) ids/geography/design variables,
* (ii) the outcome and regressor families -- here including every ingredient
* of the shift-share instrument (origin-state migrant stocks, predicted
* flows, truncation dummies) -- (iii) the primary controls, (iv) >50%
* missing, (v) the 0.8 correlation screen. This yields d = 85.
local excluded icpsrfip_1 icpsrfip year statename state_po countyname fips     ///
    office candidate party candidatevotes totalvotes version xcoord ycoord     ///
    decade gisjoin area_sqmi decade_1 gisjoin_1 area_sqmi_1 area_sqmii         ///
    statefip county South border North D unincorp_1860 clon10 clat10           ///
    km_grid_cel_code insamp Trump_share pct_Southerners_white                  ///
    votes_1940 votes_1948 votes_2000 votes_2016 tot_vote_pres                  ///
    elec1940 elecpop40 elec1948 elecpop48 elec2000 elecpop00 elec2016 elecpop16 ///
    Prep1900 FEpairRep1900 FEpairRep0040 FEpair1870share                       ///
    popdens_hist popdens_hist_00
local ylist
qui ds, has(type numeric)
foreach v in `r(varlist)' {
    local skip = 0
    if strpos(" `excluded' ", " `v' ")              local skip = 1
    if regexm("`v'", "^(D?pct_)?Southerners_white") local skip = 1
    if regexm("`v'", "^Southerners_white")          local skip = 1
    if regexm("`v'", "^D_pct_Southerners")          local skip = 1
    if regexm("`v'", "^(iv_|dummy_trun)")           local skip = 1
    if regexm("`v'", "_share_1900$")                local skip = 1
    if regexm("`v'", "^(Southerners_black|Northerners_white|pct_Southerners_black|pct_Northerners_white|Dpct_Northerners_white)_[0-9]+$") local skip = 1
    if regexm("`v'", "^pct_Northerners_white_hat")  local skip = 1
    if regexm("`v'", "^dummy_truncate")             local skip = 1
    if regexm("`v'", "^missing_")                   local skip = 1
    if strpos(" `ctrl' ", " `v' ")                  local skip = 1
    if `skip' continue
    qui count if insamp
    local nin = r(N)
    qui count if missing(`v') & insamp
    if r(N) >= 0.5*`nin' continue
    local maxc = 0
    foreach r in Trump_share pct_Southerners_white `ctrl' {
        cap qui corr `v' `r' if insamp
        if _rc continue
        if abs(r(rho)) > `maxc' local maxc = abs(r(rho))
    }
    if `maxc' < 0.8 local ylist `ylist' `v'
}
di as text "auxiliary outcomes selected: " `: word count `ylist''

tempfile bzbase
save `bzbase'

*--- (6a) TMO augmenting the original grid clustering, with diagnostics
sjlog using "$PPR/examples/bazziExample.tex", replace
tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5) ///
    plothist plothistnbins(100) plotse file("figures/bazzi")
sjlog close, replace
scalar se_bz_tmo = e(tmo_se)
scalar thres_bz  = e(threshold)
di as text "TMO/cluster ratio: " %6.3f se_bz_tmo/se_bz_base
cap erase "$PPR/figures/bazzi_dyad.dta"

*--- (6b) method comparison (feeds the paper's comparison table)
* TMO without augmenting the original clustering
use `bzbase', clear
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', robust a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5)
di as text "[TMO alone]   SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_bz_base ///
    " pct=" %6.3f e(pct_ge_thres)

* TMO augmenting state-level clusters
use `bzbase', clear
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(statefip) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5)
di as text "[TMO+state]   SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_bz_base ///
    " pct=" %6.3f e(pct_ge_thres)

* Conley 150mi only
use `bzbase', clear
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) distthreshold(150) miles thresholdoff
di as text "[Conley150]   SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_bz_base

* TMO + Conley 150mi; save the pair-level file for the predictors table
use `bzbase', clear
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) distthreshold(150) miles savedyad file("$TMP/bazzi")
di as text "[TMO+Con150]  SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_bz_base ///
    " pct=" %6.3f e(pct_ge_thres)

* SCPC and TMO + SCPC: scpc supports ivregress 2sls, so the specification is
* rewritten with explicit state dummies (identical point estimates)
use `bzbase', clear
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`ylist') i(fips) misslimit(0.5) ///
    lat(clat10) lon(clon10) ///
    scpc_cmd(ivregress 2sls Trump_share `ctrl' i.statefip ///
    dummy_trunmig_hat1_00_2 (pct_Southerners_white = iv_mig_hat1_00_2), robust)
di as text "[SCPC alone]  SE=" %7.4f scalar(scpc_se) " ratio=" %6.3f scalar(scpc_se)/se_bz_base
di as text "[TMO+SCPC]    SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_bz_base

*--- (6c) predictors of highly-correlated county pairs (feeds the paper's
*    Table-1-style diagnostic table)
use `bzbase', clear
keep fips countyname statefip km_grid_cel_code totpop
qui gduplicates drop
rename (fips countyname statefip km_grid_cel_code totpop) (id1 nm1 st1 gr1 pop1)
tempfile bz_c1
save `bz_c1'
rename (id1 nm1 st1 gr1 pop1) (id2 nm2 st2 gr2 pop2)
tempfile bz_c2
save `bz_c2'

use "$TMP/bazzi_dyad.dta", clear
keep if id1!=id2
gen byte above = (abs(corr)>=thres_bz) & !missing(corr)
merge m:1 id1 using `bz_c1', keep(3) nogen
merge m:1 id2 using `bz_c2', keep(3) nogen
gen byte near150  = dist<=150
gen byte samegrid = gr1==gr2 & !missing(gr1)
gen byte samest   = st1==st2 & !missing(st1)
gen dpop = abs(pop1-pop2)
qui _pctile dpop, p(10)
gen byte nearpop = dpop<=r(r1) if !missing(dpop)
gen byte anyclose = near150|samegrid|samest|nearpop
di as text _n "Predictors of selected county pairs (vs all pairs):"
foreach v in near150 samegrid samest nearpop anyclose {
    qui sum `v' if above
    local pa = 100*r(mean)
    qui sum `v'
    di as text "  `v': " %5.1f `pa' "%  (all pairs: " %5.1f 100*r(mean) "%)"
}
gsort -above -dist
di as text _n "Most distant highly-correlated pairs:"
li nm1 nm2 corr dist in 1/5, noobs clean

*--- (6d) sensitivity of the adjustment to the auxiliary collection.
* Row 2 trims the collection to d=60 by also excluding the secondary control
* sets of the paper's own Table 2 (1900 controls, sorting confounds, terrain)
* and the raw vote counts: an over-pruned collection whose estimated df fall
* below the recommended floor of 20.
* Row 3 is a DELIBERATE misuse example: variables measuring the treatment
* (the regressor family) must be excluded by design; diagnostics won't flag them.
use `bzbase', clear
local ctrl1900 lnpopdens_hist_00 pct_mfgempl_00 pct_btot_00 pct_popmexico_00 ///
    pct_popgerman_00 pct_popcanada_00 pct_popireland_00 pct_popitaly_00 ///
    pct_farmacres_00 pct_farmvalue_00
local confounds share_breckinridge pct_bryan_96 oil_1900 oil_1940 AnyMines ///
    cottonmed potential_ag_prod d_coa d_riv elev_mean tri_ave d_lak ///
    tye_tfe890_500k_100_l6
local yTrim
foreach v of local ylist {
    local skip = 0
    if strpos(" `ctrl1900' ", " `v' ")  local skip = 1
    if strpos(" `confounds' ", " `v' ") local skip = 1
    if regexm("`v'", "^votes_")         local skip = 1
    if !`skip' local yTrim `yTrim' `v'
}
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`yTrim') i(fips) misslimit(0.5)
di as text "[trim]    D=" e(N_outcomes) " df=" %5.1f e(dof) " thr=" %5.3f e(threshold) ///
    " ratio=" %5.3f e(tmo_se)/se_bz_base

use `bzbase', clear
local yContam `ylist' pct_Southerners_white1900 pct_Southerners_white_brdr ///
    Dpct_Southerners_white D_pct_Southerners_white_00_40
qui tmo, cmd(ivreghdfe Trump_share (pct_Southerners_white = iv_mig_hat1_00_2) ///
    `ctrl', cluster(km_grid_cel_code) a(statefip dummy_trunmig_hat1_00_2)) ///
    x(pct_Southerners_white) ylist(`yContam') i(fips) misslimit(0.5)
di as text "[contam]  D=" e(N_outcomes) " df=" %5.1f e(dof) " thr=" %5.3f e(threshold) ///
    " ratio=" %5.3f e(tmo_se)/se_bz_base

*-------------------------------------------------------------------------------
*--- (7) Step-by-step guide application: Bernini et al. (2023)
*    Fully self-contained: all 60 auxiliary outcomes come from the paper's own
*    replication package (companion paper, Appendix E.2). Published targets:
*    coef 0.10, orig SE 0.04, d=60, df=25.8, delta*=0.54, 0.70% cross-cluster
*    pairs, TMO ratio 1.37; Conley150 1.40 [9.0%], TMO+Conley 1.51 [9.6%].
*-------------------------------------------------------------------------------
use "$SJ/data/replication/dataset_wide_1.dta", clear

* county centroids (the variable -county- holds the 5-digit FIPS as string)
preserve
use "$SJ/data/maps/cb_2018_us_county_20m.dta", clear
destring GEOID, replace
rename GEOID fips
keep fips _CX _CY
tempfile cent
save `cent'
restore
gen long fips = real(county)
merge 1:1 fips using `cent', keep(1 3) nogen

* Step 1: baseline -- Table 2 Column 4 of Bernini et al. (preferred spec);
* factor-variable syntax replaces the original xi: prefix (identical estimates)
qui reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc ///
    c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc ///
    c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc ///
    c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc ///
    ibn.STATE, nocon robust cluster(judicial_divisions_id)
di as text "Bernini Table 2 Col 4: theta = " %6.4f _b[black_share60_lit_nc] ///
    "  cluster SE = " %6.4f _se[black_share60_lit_nc] "  N = " e(N)
scalar se_bern_base = _se[black_share60_lit_nc]
gen byte insamp = e(sample)

* Step 2: auxiliary outcomes, following the companion paper's E.2 rules:
* all feasible package variables, excluding (i) ids/geography/design vars,
* (ii) the outcome family, regressor family and their interactions, (iii) the
* primary controls, (iv) >50% missing, (v) |corr|>0.8 with the outcome,
* regressor, or any primary control. This yields exactly d = 60.
local excluded county countycode FIPSTATE STATE geo judicial_divisions_id ///
    literacy_nc SMD MIXED AL insamp fips _CX _CY
local refs ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 unemp60 ///
    family_less_3000 pop60 school_low urbanB60 cotton_suitability ///
    cotton_share_land1964 pro_black_county anti_black_county rep_share_1964
local ylist
local dropped
qui ds, has(type numeric)
foreach v in `r(varlist)' {
    local skip = 0
    if strpos(" `excluded' ", " `v' ") local skip = 1
    if regexm("`v'", "^(ch_)?ShareBl_") local skip = 1
    if regexm("`v'", "^diffshareblack") local skip = 1
    if regexm("`v'", "^black_share") local skip = 1
    if regexm("`v'", "_lit(_|$)") local skip = 1
    if regexm("`v'", "^dist_lit") local skip = 1
    if strpos(" `refs' ", " `v' ") local skip = 1
    if "`v'"=="ln_rep_share_1964" local skip = 1
    if `skip'==0 {
        qui sum insamp, meanonly
        local nin = r(sum)
        qui count if missing(`v') & insamp
        if r(N) < 0.5*`nin' {
            local maxc = 0
            foreach r of local refs {
                qui corr `v' `r' if insamp
                if abs(r(rho)) > `maxc' local maxc = abs(r(rho))
            }
            if `maxc' < 0.8 local ylist `ylist' `v'
            else local dropped `dropped' `v'
        }
    }
}
local d : word count `ylist'
di as text "auxiliary outcomes selected: `d'  (correlation screen removed: `dropped')"

tempfile bernbase
save `bernbase'

* Step 3: run tmo augmenting the original judicial-division clustering
sjlog using "$PPR/examples/berniniExample.tex", replace
tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc ///
    c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc ///
    c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc ///
    c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc ///
    ibn.STATE, nocon robust cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist(`ylist') i(fips) misslimit(0.5) ///
    plothist plothistnbins(100) plotse file("figures/bernini")
sjlog close, replace
scalar thres_bern = e(threshold)
di as text "TMO/cluster ratio: " %6.3f e(tmo_se)/se_bern_base
cap erase "$PPR/figures/bernini_dyad.dta"

* Steps 5-6: distance-based runs (Conley 150mi as in the companion paper) and
* the pair-level file used for the predictors table
use `bernbase', clear
qui tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc ///
    c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc ///
    c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc ///
    c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc ///
    ibn.STATE, nocon robust cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist(`ylist') i(fips) misslimit(0.5) ///
    lat(_CY) lon(_CX) distthreshold(150) miles savedyad file("$TMP/bern")
di as text "[TMO+Conley150] ratio=" %6.3f e(tmo_se)/se_bern_base ///
    " pct=" %6.3f e(pct_ge_thres)

use `bernbase', clear
qui tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
    c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc ///
    c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc ///
    c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc ///
    c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc ///
    ibn.STATE, nocon robust cluster(judicial_divisions_id)) ///
    x(black_share60_lit_nc) ylist(`ylist') i(fips) misslimit(0.5) ///
    lat(_CY) lon(_CX) distthreshold(150) miles thresholdoff
di as text "[Conley150 only] ratio=" %6.3f e(tmo_se)/se_bern_base

* Step 5: predictors of the selected pairs (feeds the paper's table)
preserve
use `bernbase', clear
keep fips STATE judicial_divisions_id pop60 family_less_3000 urbanB60
qui gduplicates drop
rename (fips STATE judicial_divisions_id pop60 family_less_3000 urbanB60) ///
       (id1 st1 jd1 pop1 pov1 urb1)
tempfile bc1
save `bc1'
rename (id1 st1 jd1 pop1 pov1 urb1) (id2 st2 jd2 pop2 pov2 urb2)
tempfile bc2
save `bc2'
restore

use "$TMP/bern_dyad.dta", clear
keep if id1!=id2
gen byte above = (abs(corr)>=thres_bern) & !missing(corr)
merge m:1 id1 using `bc1', keep(3) nogen
merge m:1 id2 using `bc2', keep(3) nogen
gen byte near150 = dist<=150
gen byte samest  = st1==st2
gen byte samejd  = jd1==jd2
foreach p in pop pov urb {
    gen d`p' = abs(`p'1-`p'2)
    qui _pctile d`p', p(10)
    gen byte near`p' = d`p'<=r(r1) if !missing(d`p')
}
gen byte anyclose = near150|samest|samejd|nearpop|nearpov|nearurb
di as text _n "Predictors of selected county pairs (vs all pairs):"
foreach v in near150 samest samejd nearpop nearpov nearurb anyclose {
    qui sum `v' if above
    local pa = 100*r(mean)
    qui sum `v'
    di as text "  `v': " %5.1f `pa' "%  (all pairs: " %5.1f 100*r(mean) "%)"
}

* Step 7: sensitivity of the adjustment to the outcome collection
use `bernbase', clear
local yfam `ylist' pop50 urban50 family_less50_2000 school_low50 ///
    familypovchange cotton_share_land1945 all_officials_1980
local ycontam `ylist' black_share40 black_share50 diffshareblack60_50
foreach variant in yfam ycontam {
    use `bernbase', clear
    qui tmo, cmd(reg ch_ShareBl_AllOfficials black_share60_lit_nc black_share60 ///
        c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc ///
        c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc ///
        c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc ///
        c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc ///
        ibn.STATE, nocon robust cluster(judicial_divisions_id)) ///
        x(black_share60_lit_nc) ylist(``variant'') i(fips) misslimit(0.5)
    di as text "[`variant'] D=" e(N_outcomes) " df=" %5.1f e(dof) ///
        " thr=" %5.3f e(threshold) " ratio=" %5.3f e(tmo_se)/se_bern_base
}
