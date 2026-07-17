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
*--- (6) Real-data application: Acemoglu et al. (2019), democracy and growth
*    A diagnostic guide for WHEN to use tmo. Replicates Table 2 Column 3
*    (preferred spec) and runs the diagnostic battery of DellaVigna et al.:
*    correlation histogram + null fit, SE across thresholds, predictors of
*    highly-correlated country pairs, and the method-comparison table.
*    Practical notes:
*      - lags are pre-computed as variables (ts-operators inside cmd() break
*        once tmo subsets to the estimation sample of an unbalanced panel)
*      - the regressor of interest (dem) is listed immediately after the
*        dependent variable in cmd()
*-------------------------------------------------------------------------------
use "$ROOT/supporting_material/replication_files_ddcg/DDCGdata_final.dta", clear
xtset wbcode2 year
forval j = 1/4 {
    gen ly`j' = l`j'.y
}

* auxiliary outcomes: country-year economic variables from the replication
* package (education, trade, investment, TFP, mortality, taxes, reforms,
* population, external assets); democracy measures are excluded
local ylist lp_bl ls_bl lh_bl taxratio secenr prienr tradewb mortnew ginv ///
    rtfpna unrest unrestn marketref gfa nfa totalassets totalliabilities ///
    nfagdp loginvpc ltfp ltrade2 lprienr lsecenr lgov lmort ///
    PopulationtotalSPPOPTOTL Populationages014oftotal ///
    Populationages1564oftota gdppercapitaconstant2000us rgdpl2 rgdpna_full

tempfile aceb
save `aceb'

* baseline: same estimates as ANRR's xtreg, fe (Table 2 Column 3)
qui areg y dem ly1 ly2 ly3 ly4 yy*, absorb(wbcode2) cluster(wbcode2)
scalar se_ace_base = _se[dem]
di as text "ANRR Table 2 Col 3:  beta = " %6.4f _b[dem] "   cluster SE = " %6.4f _se[dem]

*--- (6a) TMO augmenting the original clustering, with diagnostics
sjlog using "$PPR/examples/acemogluExample.tex", replace
tmo, cmd(areg y dem ly1 ly2 ly3 ly4 yy*, absorb(wbcode2) cluster(wbcode2)) ///
    x(dem) ylist(`ylist') i(wbcode2) t(year) ///
    plothist plothistnbins(100) plotse savedyad file("figures/acemoglu")
sjlog close, replace

scalar se_ace_tmo = e(tmo_se)
scalar thr_ace    = e(threshold)
di as text "TMO/cluster ratio: " %6.3f se_ace_tmo/se_ace_base

*--- (6b) method comparison (feeds the paper's comparison table)
use `aceb', clear
qui tmo, cmd(areg y dem ly1 ly2 ly3 ly4 yy*, absorb(wbcode2) cluster(wbcode2)) ///
    x(dem) ylist(`ylist') i(wbcode2) t(year) ///
    lat(cen_lat) lon(cen_lon) distthreshold(650) miles thresholdoff
di as text "[Conley650]  SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_ace_base

use `aceb', clear
qui tmo, cmd(areg y dem ly1 ly2 ly3 ly4 yy*, absorb(wbcode2) cluster(wbcode2)) ///
    x(dem) ylist(`ylist') i(wbcode2) t(year) ///
    lat(cen_lat) lon(cen_lon) distthreshold(650) miles
di as text "[TMO+Con650] SE=" %7.4f e(tmo_se) " ratio=" %6.3f e(tmo_se)/se_ace_base ///
    " pct=" %6.3f e(pct_ge_thres)

*--- (6c) predictors of highly-correlated country pairs (feeds the paper's
*    Table-1-style diagnostic table)
use `aceb', clear
keep wbcode2 country_name cen_lat cen_lon region gdp1960
qui gduplicates drop
rename (wbcode2 country_name cen_lat cen_lon region gdp1960) ///
       (id1 name1 lat1 lon1 reg1 gdp1)
tempfile ace_c1
save `ace_c1'
rename (id1 name1 lat1 lon1 reg1 gdp1) (id2 name2 lat2 lon2 reg2 gdp2)
tempfile ace_c2
save `ace_c2'

use "$PPR/figures/acemoglu_dyad.dta", clear
keep if id1!=id2
merge m:1 id1 using `ace_c1', keep(3) nogen
merge m:1 id2 using `ace_c2', keep(3) nogen
qui geodist lat1 lon1 lat2 lon2, gen(distmi) miles sphere
gen byte above    = (abs(corr)>=thr_ace) & !missing(corr)
gen byte near650  = distmi<=650
gen byte samereg  = (reg1==reg2) & !missing(reg1)
gen dgdp = abs(gdp1-gdp2)
qui _pctile dgdp, p(10)
gen byte neargdp  = dgdp<=r(r1) if !missing(dgdp)
gen byte anyclose = near650==1 | samereg==1 | neargdp==1

di as text _n "Predictors of |corr|>=threshold country pairs (vs all pairs):"
foreach v in near650 samereg neargdp anyclose {
    qui sum `v' if above==1
    local pa = 100*r(mean)
    qui sum `v'
    di as text "  `v': " %5.1f `pa' "%  (all pairs: " %5.1f 100*r(mean) "%)"
}
gsort -above -distmi
di as text _n "Most distant highly-correlated pairs:"
li name1 name2 corr distmi in 1/5, noobs clean
erase "$PPR/figures/acemoglu_dyad.dta"

*--- (6d) sensitivity of the adjustment to the auxiliary collection (feeds
*    the paper's sensitivity table). Row 3 is a DELIBERATE misuse example:
*    treatment measures must be excluded by design, diagnostics won't flag them
use `aceb', clear
local yThin lp_bl ls_bl lh_bl taxratio tradewb mortnew ginv rtfpna unrest ///
    marketref nfagdp lgov prienr secenr PopulationtotalSPPOPTOTL ///
    Populationages014oftotal Populationages1564oftota
qui tmo, cmd(areg y dem ly1 ly2 ly3 ly4 yy*, absorb(wbcode2) cluster(wbcode2)) ///
    x(dem) ylist(`yThin') i(wbcode2) t(year)
di as text "[thin]    D=" e(N_outcomes) " df=" %5.1f e(dof) " thr=" %5.3f e(threshold) ///
    " ratio=" %5.3f e(tmo_se)/se_ace_base

use `aceb', clear
local yContam `ylist' demFH demPOL demBMR demCGV polity2 demevent revevent
qui tmo, cmd(areg y dem ly1 ly2 ly3 ly4 yy*, absorb(wbcode2) cluster(wbcode2)) ///
    x(dem) ylist(`yContam') i(wbcode2) t(year)
di as text "[contam]  D=" e(N_outcomes) " df=" %5.1f e(dof) " thr=" %5.3f e(threshold) ///
    " ratio=" %5.3f e(tmo_se)/se_ace_base

*-------------------------------------------------------------------------------
*--- (7) Step-by-step guide application: Bernini et al. (2023)
*    Fully self-contained: all 60 auxiliary outcomes come from the paper's own
*    replication package (companion paper, Appendix E.2). Published targets:
*    coef 0.10, orig SE 0.04, d=60, df=25.8, delta*=0.54, 0.70% cross-cluster
*    pairs, TMO ratio 1.37; Conley150 1.40 [9.0%], TMO+Conley 1.51 [9.6%].
*-------------------------------------------------------------------------------
use "$ROOT/supporting_material/Bernini et al (2023)/datasets/dataset_wide_1.dta", clear

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
