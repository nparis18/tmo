/*------------------------------------------------------------------------------

Paper:			- Race, Representation and Local Governments in the U.S. South:
					The Effect of the Voting Rights Act

Journal:		- Journal of Political Economy, 2022

Authors:		- Andrea Bernini:
					University of Oxford
				- Giovanni Facchini:
					University of Nottingham, CEPR, CES-Ifo, GEP and IZA
				- Cecilia Testa:
					University of Nottingham, LdA and NICEP

Instructions:	- README.pdf

Input files:	- countiescoord.dta
				- dataset_wide_1.dta
				- dataset_wide_2.dta
				- dataset_wide_3.dta
				- dataset_wide_4.dta
				- dataset_wide_5.dta
				- dataset_panel_1.dta
				- dataset_panel_2.dta
				- dataset_panel_3.dta
				- dataset_border_1.dta
				- dataset_border_2.dta

Output files:   - Figures 1 to 10
				- Tables 1 to 7

------------------------------------------------------------------------------*/


clear all

// Install the graphics scheme Tufte
ssc install scheme_tufte
set scheme tufte

// Change the directory
global USER					"C:\Users\main_folder"

global dataset				"$USER\data_code\datasets"
global build				"$USER\build_figures"
global save_figures_tables	"$USER\ms"

// .ado files
run "$USER\data_code\bernini_facchini_testa_2022_twowaycluster_1.ado"
run "$USER\data_code\bernini_facchini_testa_2022_twowaycluster_2.ado"

// Global
global controls_event "urbanB60 unemp60 family_less_3000 pop60 school_low cotton_suitability pro_black_county anti_black_county"
global event "black_share60_lit_nc_d1964 black_share60_lit_nc_d1968 black_share60_lit_nc_d1970 black_share60_lit_nc_d1972 black_share60_lit_nc_d1974 black_share60_lit_nc_d1976 black_share60_lit_nc_d1978 black_share60_lit_nc_d1980"
global event_2 "black_share60_lit_nc_d1960 black_share60_lit_nc_d1964 black_share60_lit_nc_d1968 black_share60_lit_nc_d1972 black_share60_lit_nc_d1976 black_share60_lit_nc_d1980"
global event_3 "black_share60_lit_nc_d1962 black_share60_lit_nc_d1967 black_share60_lit_nc_d1972 black_share60_lit_nc_d1977 black_share60_lit_nc_d1982"
global control_pre "urbanB60 unemp60 family_less_3000 pop60 school_low cotton_suitability"
global pre_cr_trend "ch_kkk_66_40_pcw ch_black_lynching60_40 ch_cotton_sh_land_64_45 ch_naacp_64_42_pcb lndiff_rep_share_64_40 lndiff_rep_share_60_40 lndiff_pop_50_60 ch_black_share_60_50"
global pre_cr_level_i "kkk_klavern_15_40_pcw black_lynching_1930_1940 cotton_share_land1945 naacp_branch_1942_pcb ln_rep_share_1952 ln_rep_share_1940"
global pre_cr_level_f "kkk_klavern_64_66_pcw black_lynching_1950_1964 cotton_share_land1964 naacp_branch_1964_pcb ln_rep_share_1964 ln_rep_share_1960"
global pre_pol_trend "lndiff_tnt_pres_60_40 lndiff_tnt_gov_60_40 lndiff_winner_gov lndiff_lower_60_50 lndiff_upper_60_50"
global pre_pol_level_i "ln_tnt_pres_tot_1940 ln_tnt_gov_1940 ln_winner_gov_40 ln_lower_50 ln_upper_50"
global pre_pol_level_f "ln_tnt_pres_tot_1960 ln_tnt_gov_1960 ln_winner_gov_60 ln_lower_60 ln_upper_60"
global interaction_full "c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc"
global interaction_full_SMD "c.urbanB60#literacy_nc#SMD c.unemp60#literacy_nc#SMD c.family_less_3000#literacy_nc#SMD c.pop60#literacy_nc#SMD c.school_low#literacy_nc#SMD c.cotton_suitability#literacy_nc#SMD c.cotton_share_land1964#literacy_nc#SMD c.anti_black_county#literacy_nc#SMD c.pro_black_county#literacy_nc#SMD c.rep_share_1964#literacy_nc#SMD"



********************************************************************************
* FIGURE 1
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

spmap all_officials_1964 using "countiescoord.dta", ///
id(geo) ndfcolor(gs14) legstyle(2) legjunction(" - ") ///
legend(ring(1) pos(6) size(*1.8)) title("") ndlabel("No data") ///
fcolor(white gs9 gs5 gs2) clmethod(custom) ///
clbreaks(0 0.9 1.9 2.9 20) legorder(lohi) legend(label(2 "0") label(3 "1") ///
label(4 "2") label(5 "3+") color(dknavy*2.0))

cd $save_figures_tables
graph export "figure_1.eps", replace



********************************************************************************
* FIGURE 2
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

spmap all_officials_1980 using "countiescoord.dta", ///
id(geo) ndfcolor(gs14) legstyle(2) legjunction(" - ") ///
legend(ring(1) pos(6) size(*1.8)) title("") ndlabel("No data") ///
fcolor(white gs9 gs5 gs2) clmethod(custom) ///
clbreaks(0 0.9 1.9 2.9 56) legorder(lohi) legend(label(2 "0") label(3 "1") ///
label(4 "2") label(5 "3+") color(dknavy*2.0))	

cd $save_figures_tables
graph export "figure_2.eps", replace



********************************************************************************
* FIGURE 3
********************************************************************************

// All officials
cd $dataset
use "dataset_panel_1.dta", clear
		
drop if ShareBl_AllOfficials == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

foreach x in 62 64 68 70 72 74 76 78 80 {

xi, prefix(Tr): reg ShareBl_AllOfficials black_share60 ///
c.black_share60#literacy_nc $controls_event ///
i.FIPSTATE if year==19`x' , ///
robust cluster (judicial_divisions_id)
margins, over(literacy_nc) dydx(black_share60) post
est sto sharebl`x'

}

coefplot (sharebl62, label(1962) msymbol(O) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl64, label(1964) msymbol(Oh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl68, label(1968) msymbol(S) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl70, label(1970) msymbol(Sh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl72, label(1972) msymbol(D) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl74, label(1974) msymbol(Dh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl76, label(1976) msymbol(T) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl78, label(1978) msymbol(Th) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl80, label(1980) msymbol(X) msize(medlarge) ///
color(black) ciopts(lcolor(black))) , ///
subtitle(, fcolor(white)) xlabel("") yscale(range(0 0.4)) ///
ylabel(0(0.1)0.4, format(%9.1f)) drop($controls_event TrFIPSTATE_*_cons) ///
vertical yline(0,  lpattern(dash)) bycoefs legend(rows(2) ) ///
byopts(title("All Local Officials")) ylabel(,labsize(medium)) ///
graphregion(color(white)) bgcol(color(white)) omitted

cd $build
graph save "figure_3_a.gph", replace

// County governments
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllCountyOfficials == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

foreach x in 62 64 68 70 72 74 76 78 80 {

xi, prefix(Tr): reg ShareBl_AllCountyOfficials black_share60 ///
c.black_share60#literacy_nc $controls_event ///
i.FIPSTATE if year==19`x' , ///
robust cluster (judicial_divisions_id)
margins, over(literacy_nc) dydx(black_share60) post
est sto sharebl`x'

}

coefplot (sharebl62, label(1962) msymbol(O) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl64, label(1964) msymbol(Oh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl68, label(1968) msymbol(S) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl70, label(1970) msymbol(Sh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl72, label(1972) msymbol(D) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl74, label(1974) msymbol(Dh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl76, label(1976) msymbol(T) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl78, label(1978) msymbol(Th) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl80, label(1980) msymbol(X) msize(medlarge) ///
color(black) ciopts(lcolor(black))) , ///
subtitle(, fcolor(white)) xlabel("") yscale(range(0 0.4)) ///
ylabel(0(0.1)0.4, format(%9.1f)) drop($controls_event TrFIPSTATE_* _cons) ///
vertical yline(0, lpattern(dash)) bycoefs legend(rows(2)) ///
byopts(title("County Governments")) ylabel(,labsize(medium)) omitted

cd $build
graph save "figure_3_b.gph", replace

// School district
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllEducation == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

foreach x in 62 64 68 70 72 74 76 78 80 {

xi, prefix(Tr): reg ShareBl_AllEducation black_share60 ///
c.black_share60#literacy_nc $controls_event ///
i.FIPSTATE if year==19`x' , ///
robust cluster (judicial_divisions_id)
margins, over(literacy_nc) dydx(black_share60) post
est sto sharebl`x'

}

coefplot (sharebl62, label(1962) msymbol(O) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl64, label(1964) msymbol(Oh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl68, label(1968) msymbol(S) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl70, label(1970) msymbol(Sh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl72, label(1972) msymbol(D) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl74, label(1974) msymbol(Dh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl76, label(1976) msymbol(T) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl78, label(1978) msymbol(Th) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl80, label(1980) msymbol(X) msize(medlarge) ///
color(black) ciopts(lcolor(black))) , ///
subtitle(, fcolor(white)) xlabel("") yscale(range(0 0.4)) ///
ylabel(0(0.1)0.4, format(%9.1f)) drop($controls_event TrFIPSTATE_* _cons) ///
vertical yline(0, lpattern(dash)) bycoefs legend(rows(2)) ///
byopts(title("School Districts")) ylabel(,labsize(medium)) omitted

cd $build
graph save "figure_3_c.gph", replace

// Municipality
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllMunicipality == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

foreach x in 62 64 68 70 72 74 76 78 80 {

xi, prefix(Tr): reg ShareBl_AllMunicipality black_share60 ///
c.black_share60#literacy_nc $controls_event ///
i.FIPSTATE if year==19`x' , ///
robust cluster (judicial_divisions_id)
margins, over(literacy_nc) dydx(black_share60) post
est sto sharebl`x'

}

coefplot (sharebl62, label(1962) msymbol(O) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl64, label(1964) msymbol(Oh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl68, label(1968) msymbol(S) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl70, label(1970) msymbol(Sh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl72, label(1972) msymbol(D) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl74, label(1974) msymbol(Dh) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl76, label(1976) msymbol(T) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl78, label(1978) msymbol(Th) msize(medlarge) ///
color(black) ciopts(lcolor(black))) ///
(sharebl80, label(1980) msymbol(X) msize(medlarge) ///
color(black) ciopts(lcolor(black))) , ///
subtitle(, fcolor(white)) xlabel("") yscale(range(0 0.4)) ///
ylabel(0(0.1)0.4, format(%9.1f)) drop($controls_event TrFIPSTATE_* _cons) ///
vertical yline(0, lpattern(dash)) bycoefs legend(rows(2)) ///
byopts(title("Municipalities")) ylabel(,labsize(medium)) omitted

cd $build
graph save "figure_3_d.gph", replace

// Combine all the plots
cd $build
graph combine "figure_3_a.gph" "figure_3_b.gph" ///
"figure_3_c.gph" "figure_3_d.gph" , ///
graphregion(color(white)) ycommon

cd $save_figures_tables
graph export "figure_3.eps", replace



********************************************************************************
* FIGURE 4
********************************************************************************

// All Officials
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllOfficials == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

xi, prefix(Tr): xtreg ShareBl_AllOfficials ///
black_share60_d1964-cotton_suitability_d1980 ///
black_share60_lit_nc_d1964-cotton_suitability_lit_nc_d1980 ///
i.FIPSTATE*i.year , ///
fe robust cluster(judicial_divisions_id)
est sto ShareBl_AllOfficials

coefplot (ShareBl_AllOfficials, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.1f)) ///
title(All Local Officials) graphregion(color(white))

cd $build
graph save "figure_4_a.gph", replace

// County Governments
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllCountyOfficials == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

xi, prefix(Tr): xtreg ShareBl_AllCountyOfficials ///
black_share60_d1964-cotton_suitability_d1980 ///
black_share60_lit_nc_d1964-cotton_suitability_lit_nc_d1980 ///
i.FIPSTATE*i.year , ///
fe robust cluster(judicial_divisions_id) 
est sto AllCountyOfficials

coefplot (AllCountyOfficials, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.1f)) ///
title(County Governments) graphregion(color(white))

cd $build
graph save "figure_4_b.gph", replace

// School District
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllEducation == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

xi, prefix(Tr): xtreg ShareBl_AllEducation ///
black_share60_d1964-cotton_suitability_d1980 ///
black_share60_lit_nc_d1964-cotton_suitability_lit_nc_d1980 ///
i.FIPSTATE*i.year , ///
fe robust cluster(judicial_divisions_id)
est sto ShareBl_AllEducation

coefplot (ShareBl_AllEducation, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.1f)) ///
title(School Districts) graphregion(color(white))

cd $build
graph save "figure_4_c.gph", replace

// Municipality
cd $dataset
use "dataset_panel_1.dta", clear

drop if ShareBl_AllMunicipality == .
bysort countycode: gen nyear = [_N]
keep if nyear == 9

xi, prefix(Tr): xtreg ShareBl_AllMunicipality ///
black_share60_d1964-cotton_suitability_d1980 ///
black_share60_lit_nc_d1964-cotton_suitability_lit_nc_d1980 ///
i.FIPSTATE*i.year , ///
fe robust cluster(judicial_divisions_id) 
est sto ShareBl_AllMunicipality

coefplot (ShareBl_AllMunicipality, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.1f)) ///
title(Municipalities) graphregion(color(white))

cd $build
graph save "figure_4_d.gph", replace

// Combine all the plots
cd $build
graph combine "figure_4_a.gph" "figure_4_b.gph" ///
"figure_4_c.gph" "figure_4_d.gph" , ///
graphregion(color(white)) ycommon imargin(1 1 1 1)

cd $save_figures_tables
graph export "figure_4.eps", replace



********************************************************************************
* FIGURE 5
********************************************************************************

cd $dataset
use "dataset_wide_3.dta", clear

spmap ch_incumbency using "countiescoord.dta", ///
id(geo) ndfcolor(white) legstyle(2) legjunction(" - ") ///
legend(ring(1) pos(6) size(*1.8)) title("") ndlabel("No data") ///
fcolor(gs14 gs11 gs8 gs5 gs3 black) ///
clmethod(custom) clbreaks(-1.1 -.5001 -.2501 .0001 .2501 .5001 1.1) ///
legorder(lohi) legend(label(2 "Below -50pp") label(3 "-50pp to -25pp") ///
label(4 "-25pp to 0pp") label(5 "0pp to 25pp") label(6 "25pp to 50pp") ///
label(7 "Above 50pp") color(dknavy*2.0))

cd $save_figures_tables
graph export "figure_5.eps", replace



********************************************************************************
* FIGURE 6
********************************************************************************

cd $dataset
use "dataset_panel_2.dta", clear

// Incumbency
xi, prefix(Tr): xtreg incumbency ///
black_share60_d1964-cotton_suitability_d1980 ///
black_share60_lit_nc_d1964-cotton_suitability_lit_nc_d1980 ///
i.FIPSTATE*i.year , ///
fe robust cluster(judicial_divisions_id)
est sto incumbency

coefplot (incumbency, msymbol(O) msize(large) color(black) ///
ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event_2) ///
xlabels( , labs() labgap()) ///
ylabels( -0.4(0.2)0.2, labs() labgap() format(%9.1f)) ///
title(Incumbency Rate, size()) graphregion(color(white))

cd $build
graph save "figure_6_a.gph", replace

// County Governments
xi, prefix(Tr): xtreg ShareBl_CountyGoverningBody ///
black_share60_d1964-cotton_suitability_d1980 ///
black_share60_lit_nc_d1964-cotton_suitability_lit_nc_d1980 ///
i.FIPSTATE*i.year if incumbency != . , ///
fe robust cluster(judicial_divisions_id)
est sto ShareBl_AllCountyOfficials

coefplot (ShareBl_AllCountyOfficials, msymbol(O) msize(large) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event_2) ///
xlabels( , labs() labgap()) ///
ylabels( 0(0.1)0.2, labs() labgap() format(%9.1f)) ///
title(Share of BEO, size()) graphregion(color(white))
	
cd $build
graph save "figure_6_b.gph", replace

// Combine all the plots
cd $build
graph combine "figure_6_a.gph" "figure_6_b.gph", ///
fysize() graphregion(fcolor(gs16)) col(1)

cd $save_figures_tables
graph export "figure_6.eps", replace



********************************************************************************
* FIGURE 7
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

drop if dist_lit_nolit_25 == . | dist_lit_nolit_50 == . | ///
dist_lit_nolit_75 == . | dist_lit_nolit_100 == .

ttest black_share60 if dist_lit_nolit_25 > 0, by(literacy_nc)

local mean1 = r(mu_1)
local mean2 = r(mu_2)

local se = r(se)
gen SEdiffmean1 = r(se)
gen diffmean1 = r(mu_1) - r(mu_2)

ttest black_share60 if dist_lit_nolit_50 > 0, by(literacy_nc)

local mean1 = r(mu_1)
local mean2 = r(mu_2)

local se = r(se)
gen SEdiffmean2 = r(se)
gen diffmean2 = r(mu_1) - r(mu_2)

ttest black_share60 if dist_lit_nolit_75 > 0, by(literacy_nc)

local mean1 = r(mu_1)
local mean2 = r(mu_2)

local se = r(se)
gen SEdiffmean3 = r(se)
gen diffmean3 = r(mu_1) - r(mu_2)

ttest black_share60 if dist_lit_nolit_100 > 0, by(literacy_nc)

local mean1 = r(mu_1)
local mean2 = r(mu_2)

local se = r(se)
gen SEdiffmean4 = r(se)
gen diffmean4 = r(mu_1) - r(mu_2)

foreach x in 25 50 75 100 {
gen dist`x' = `x'
}

collapse dist25 dist50 dist75 dist100 dist_lit_nolit_* SEdiffmean* diffmean*

stack dist25 dist_lit_nolit_25 SEdiffmean1 diffmean1 ///
dist50 dist_lit_nolit_50 SEdiffmean2 diffmean2 ///
dist75 dist_lit_nolit_75 SEdiffmean3 diffmean3 ///
dist100 dist_lit_nolit_100 SEdiffmean4 diffmean4 ///
, into (distance avgdist sediff_means diff_means) clear

serrbar diff_means sediff_means distance, ///
scale(2) mvopts(mcolor(black) msize(medium) msymbol(circle) mlabel() ///
mlabposition(2) mlabcolor(black) mlabformat(%3.2f)) ///
ytitle(Difference in Means) yscale(range(-15 7)) ylabel(-15 (5) 7) ///
xtitle(Distance from the Border) xscale(range(15 110)) xlabel(25(25)100) ///
yline(0, lpattern(dash)) lwidth(small) lcolor(black) ///
graphregion(color(white)) msize(*0)

cd $save_figures_tables
graph export "figure_7.eps", replace



********************************************************************************
* FIGURE 8
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

foreach var of varlist black_share60 black_share40 black_share50 pop60 pop50 ///
unemp60 unemp50 family_less_3000 family_less50_2000 urbanB60 urban50 ///
school_low school_low50 cotton_suitability pro_black_county anti_black_county {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc , robust cluster(judicial_divisions_id)
est sto `var'_std

}

coefplot (black_share60_std, label(Black Share) msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(unemp60_std, label(Unemployment) msymbol(Oh) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(family_less_3000_std, label(Poverty) msymbol(S) msize(medium) ///
color(black) ciopts(lcolor(black)) ) ///
(pop60_std, label(Population) msymbol(Sh) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(urbanB60_std, label(Percent Urban) msymbol(D) msize(medium) ///
color(black) ciopts(lcolor(black)))  ///
(school_low_std, label(Unskilled) msymbol(Dh) msize(medium) ///
color(black) ciopts(lcolor(black))) /// 
(cotton_suitability_std, label(Cotton Suitability) msymbol(T) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(pro_black_county_std, label(Pro-black) msymbol(Th) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(anti_black_county_std, label(Anti-black) msymbol(X) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc) title(Overall Sample) levels(95) ytitle("") ///
ylabel("") legend(rows(3) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))

cd $build
graph save "figure_8_a.gph", replace

// Border dataset
cd $dataset
use "dataset_border_1.dta", clear

foreach var of varlist black_share60 black_share40 black_share50 pop60 pop50 ///
unemp60 unemp50 family_less_3000 family_less50_2000 urbanB60 urban50 ///
school_low school_low50 cotton_suitability pro_black_county anti_black_county {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc i.pair_id, robust cluster(judicial_border)
est sto `var'_std

}

coefplot (black_share60_std, label(Black Share) msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(unemp60_std, label(Unemployment) msymbol(Oh) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(family_less_3000_std, label(Poverty) msymbol(S) msize(medium) ///
color(black) ciopts(lcolor(black)) ) ///
(pop60_std, label(Population) msymbol(Sh) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(urbanB60_std, label(Percent Urban) msymbol(D) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(school_low_std, label(Unskilled) msymbol(Dh) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(cotton_suitability_std, label(Cotton Suitability) msymbol(T) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(pro_black_county_std, label(Pro-black) msymbol(Th) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
(anti_black_county_std, label(Anti-black) msymbol(X) msize(medium) ///
color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc) title(Border Sample) levels(95) ytitle("") ///
ylabel("") legend(rows(3)size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))

cd $build
graph save "figure_8_b.gph", replace

// Combine all the plots
cd $build
graph combine "figure_8_a.gph" "figure_8_b.gph", ///
graphregion(color(white)) xcommon

cd $save_figures_tables
graph export "figure_8.eps", replace



********************************************************************************
* FIGURE 9
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

// Pre-VRA trend
foreach var of varlist unempchange60_50 ruralchange60_50 urbanchange60_50 ///
familypovchange60_50 diffshareblack60_50 lndiff_pop_50_60 ///
school_lowchange60_50 {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc , robust cluster(judicial_divisions_id)
est sto `var'_std

}

coefplot (diffshareblack60_50_std, label({&Delta} Black share) ///
msymbol(O) msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange60_50_std, label({&Delta} Unemployment) ///
msymbol(Oh) msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange60_50_std, label({&Delta} Poverty) ///
msymbol(S) msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange60_50_std, label({&Delta} Unskilled) ///
msymbol(Sh) msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange60_50_std, label({&Delta} Percent Urban) ///
msymbol(D) msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_50_60_std, label({&Delta} Population) ///
msymbol(Dh) msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc) title(Pre-VRA: Overall) levels(95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))

cd $build
graph save "figure_9_a.gph", replace

// Post-VRA trend
foreach var of varlist unempchange ruralchange urbanchange familypovchange ///
diffshareblack80 lndiff_pop_60_80 school_lowchange {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc , robust cluster(judicial_divisions_id)
est sto `var'_std

}

coefplot (diffshareblack80_std, label({&Delta} Black share) ///
msymbol(O) msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange_std, label({&Delta} Unemployment) ///
msymbol(Oh) msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange_std, label({&Delta} Poverty) ///
msymbol(S) msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange_std, label({&Delta} Unskilled) ///
msymbol(Sh) msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange_std, label({&Delta} Percent Urban) ///
msymbol(D) msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_60_80_std , label({&Delta} Population) ///
msymbol(Dh) msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc) title(Post-VRA: Overall) levels(95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))

cd $build
graph save "figure_9_b.gph", replace

// Border dataset
cd $dataset
use "dataset_border_1.dta", clear

// Pre-VRA trend
foreach var of varlist unempchange60_50 ruralchange60_50 urbanchange60_50 ///
familypovchange60_50 diffshareblack60_50 lndiff_pop_50_60 ///
school_lowchange60_50 {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc ,  robust cluster(judicial_border)
est sto `var'_std

}

coefplot (diffshareblack60_50_std, label({&Delta} Black share) ///
msymbol(O) msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange60_50_std, label({&Delta} Unemployment) ///
msymbol(Oh) msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange60_50_std, label({&Delta} Poverty) ///
msymbol(S) msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange60_50_std, label({&Delta} Unskilled) ///
msymbol(Sh) msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange60_50_std, label({&Delta} Percent Urban) ///
msymbol(D) msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_50_60_std, label({&Delta} Population) ///
msymbol(Dh) msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc) title(Pre-VRA: Border) levels(95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))

cd $build
graph save "figure_9_c.gph", replace

// Post-VRA trend
foreach var of varlist unempchange ruralchange urbanchange familypovchange ///
diffshareblack80 lndiff_pop_60_80 school_lowchange {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc , robust cluster(judicial_border)
est sto `var'_std 

}

coefplot (diffshareblack80_std, label({&Delta} Black share) ///
msymbol(O) msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange_std, label({&Delta} Unemployment) ///
msymbol(Oh) msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange_std, label({&Delta} Poverty) ///
msymbol(S) msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange_std, label({&Delta} Unskilled) ///
msymbol(Sh) msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange_std, label({&Delta} Percent Urban) ///
msymbol(D) msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_60_80_std , label({&Delta} Population) ///
msymbol(Dh) msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc) title(Post-VRA: Border) levels(95) ytitle("") ///
ylabel("") legend(rows(2)  size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))

cd $build
graph save "figure_9_d.gph", replace

// Combine all the plots
cd $build
graph combine "figure_9_a.gph" "figure_9_b.gph" ///
"figure_9_c.gph" "figure_9_d.gph", ///
graphregion(color(white)) xcommon ycommon

cd $save_figures_tables
graph export "figure_9.eps", replace



********************************************************************************
* FIGURE 10
********************************************************************************

// Current
cd $dataset
use "dataset_panel_3.dta", clear

drop if ln_county_current == .
bysort countycode: gen nyear = [_N]
keep if nyear == 6

xi, prefix(Tr): xtreg ln_county_current ///
black_share60_d1962-family_less_d1982 ///
black_share60_lit_nc_d1962-family_less_lit_nc_d1982 ///
i.FIPSTATE*i.year , ///
fe robust cluster(judicial_divisions_id)
est sto countycurr

coefplot (countycurr, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event_3) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.2f)) ///
title(Current) graphregion(color(white))

cd $build
graph save "figure_10_a.gph", replace

// Capital
cd $dataset
use "dataset_panel_3.dta", clear

drop if ln_county_cap == .
bysort countycode: gen nyear = [_N]
keep if nyear == 6

xi, prefix(Tr): xtreg ln_county_cap ///
black_share60_d1962-family_less_d1982 ///
black_share60_lit_nc_d1962-family_less_lit_nc_d1982 ///
i.FIPSTATE*i.year, ///
fe robust cluster(judicial_divisions_id)
est sto countycap

coefplot (countycap, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event_3) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.2f)) ///
title(Capital) graphregion(color(white))

cd $build
graph save "figure_10_b.gph", replace
 
// Capital, SMD
xi, prefix(Tr): xtreg ln_county_cap ///
black_share60_d1962-family_less_d1982 ///
black_share60_lit_nc_d1962-family_less_lit_nc_d1982 ///
i.FIPSTATE*i.year if SMD==1 , ///
fe robust cluster(judicial_divisions_id)
est sto countycap_SMD

coefplot (countycap_SMD, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event_3) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.2f)) ///
title(Capital SMD) graphregion(color(white))

cd $build
graph save "figure_10_c.gph", replace

// Capital, non-SMD
xi, prefix(Tr): xtreg ln_county_cap ///
black_share60_d1962-family_less_d1982 ///
black_share60_lit_nc_d1962-family_less_lit_nc_d1982 ///
i.FIPSTATE*i.year if SMD==0 , ///
fe robust cluster(judicial_divisions_id)
est sto countycap_no_SMD

coefplot (countycap_no_SMD, msymbol(O) msize(medium) ///
color(black) ciopts(lcolor(black))) , ///
vertical yline(0, lpattern(dash)) xline() keep($event_3) ///
xlabels( , labs(medium) labgap()) ///
ylabels( , labs(medium) labgap() format(%9.2f)) ///
title(Capital non-SMD) graphregion(color(white))

cd $build
graph save "figure_10_d.gph", replace

// Combine all the plots
cd $build
graph combine "figure_10_a.gph" "figure_10_b.gph" ///
"figure_10_c.gph" "figure_10_d.gph" , ///
graphregion(color(white)) imargin(2 2 2 2) ycommon

cd $save_figures_tables
graph export "figure_10.eps", replace



********************************************************************************
* TABLE 1
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

// Pre-trends
foreach var of varlist $pre_cr_trend {

xi: reg `var' black_share60_lit_nc black_share60 $control_pre ibn.STATE , ///
robust cluster(judicial_divisions_id)
eststo `var'_t		

}

foreach var of varlist  $pre_cr_level_i {

xi: reg `var' black_share60_lit_nc black_share60 $control_pre ibn.STATE , ///
robust cluster(judicial_divisions_id)
eststo `var'_i

}

foreach var of varlist  $pre_cr_level_f {

xi: reg `var' black_share60_lit_nc black_share60 $control_pre ibn.STATE , ///
robust cluster(judicial_divisions_id)
eststo `var'_f

}

// Pre-trends
foreach var of varlist $pre_pol_trend {

xi: reg `var' black_share60_lit_nc black_share60 $control_pre ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo `var'_t

}

foreach var of varlist $pre_pol_level_i {

xi: reg `var' black_share60_lit_nc black_share60 $control_pre ibn.STATE , ///
robust cluster(judicial_divisions_id)
eststo `var'_i

}

foreach var of varlist $pre_pol_level_f {

xi: reg `var' black_share60_lit_nc black_share60 $control_pre ibn.STATE , ///
robust cluster(judicial_divisions_id)
eststo `var'_f

}

// Border sample
cd $dataset
use "dataset_border_1.dta", clear

// Pre-trends
foreach var of varlist $pre_cr_trend {

twowayclusterBTF_withconst `var' ///
black_share60_lit_nc black_share60 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127, ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo `var'_b_t

}

// Pre-trends
foreach var of varlist $pre_pol_trend {

twowayclusterBTF_withconst `var' ///
black_share60_lit_nc black_share60 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127, ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo `var'_b_t

}

cd $save_figures_tables
estout ch_kkk_66_40_pcw_t ch_black_lynching60_40_t ch_cotton_sh_land_64_45_t ch_naacp_64_42_pcb_t lndiff_rep_share_64_40_t lndiff_rep_share_60_40_t lndiff_winner_gov_t lndiff_tnt_pres_60_40_t  lndiff_tnt_gov_60_40_t   lndiff_lower_60_50_t lndiff_upper_60_50_t lndiff_pop_50_60_t ch_black_share_60_50_t using table_1.tex, stats(r2_a  N, fmt(4 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n  0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(4))" " se(par label( ) fmt(4))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
title(Dependent Variable: Pre-VRA trends and levels in civil rights activism, racial composition, attitudes and political participation \label{pretrend1}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{adjustbox}{scale=0.65} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c c c c c c c c c c c c} \toprule & KKK & Lynching & Cotton & NAACP & Goldwater & Republican & Governor win & President tnt & Governor tnt & State House & State Senate & Population & Percent Black\\ & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) & (11) & (12) & (13) \\ \midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel A: Trend (Overall Sample)}} \\ \midrule")  ///
prefoot("& & & & & & & & & & & & & \\") ///
postfoot( "\midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel B: Final Level (Overall Sample)}} \\ \midrule") ///

estout kkk_klavern_64_66_pcw_f black_lynching_1950_1964_f cotton_share_land1964_f naacp_branch_1964_pcb_f ln_rep_share_1964_f ln_rep_share_1960_f ln_winner_gov_60_f ln_tnt_pres_tot_1960_f ln_tnt_gov_1960_f ln_lower_60_f ln_upper_60_f using table_1.tex, stats(r2_a  N, fmt(4 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(4))" " se(par label( ) fmt(4))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
prefoot("& & & & & & & & & & & & & \\") ///
postfoot("\midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel C: Initial Level (Overall Sample)}} \\ \midrule") ///
			
estout kkk_klavern_15_40_pcw_i black_lynching_1930_1940_i cotton_share_land1945_i naacp_branch_1942_pcb_i ln_rep_share_1952_i ln_rep_share_1940_i ln_winner_gov_40_i ln_tnt_pres_tot_1940_i  ln_tnt_gov_1940_i ln_lower_50_i ln_upper_50_i using table_1.tex, stats(r2_a  N, fmt(4 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(4))" " se(par label( ) fmt(4))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
prefoot("& & & & & & & & & & & & &  \\") ///
postfoot("\midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel D: Trend (Border Sample)}} \\ \midrule") ///

estout  ch_kkk_66_40_pcw_b_t ch_black_lynching60_40_b_t ch_cotton_sh_land_64_45_b_t ch_naacp_64_42_pcb_b_t lndiff_rep_share_64_40_b_t lndiff_rep_share_60_40_b_t lndiff_winner_gov_b_t lndiff_tnt_pres_60_40_b_t lndiff_tnt_gov_60_40_b_t lndiff_lower_60_50_b_t lndiff_upper_60_50_b_t lndiff_pop_50_60_b_t ch_black_share_60_50_b_t using table_1.tex, stats(r2_a  N, fmt(4 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(4))" " se(par label( ) fmt(4))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
prefoot("& & & & & & & & & & & & & \\") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions in Panels A-B-C and by judicial divisions and border segments in Panel D. State trends are included in Panels A-B-C and county-pair and coverage trends are included in Panel D. Controls in Panels A-B-C: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes}\end{threeparttable} \end{adjustbox} \end{center} \end{table}") 



********************************************************************************
* TABLE 2
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 , ///
nocon robust cluster(judicial_divisions_id)
eststo base1
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 ///
unemp60 family_less_3000 pop60 school_low urbanB60 cotton_suitability ///
cotton_share_land1964 ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo base2
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 ///
unemp60 pop60 family_less_3000 school_low urbanB60 cotton_suitability ///
cotton_share_land1964 pro_black_county anti_black_county rep_share_1964 ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo base3
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo base4

cd $save_figures_tables
estout base1 base2 base3 base4 ///
using table_2.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n  0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc black_share60) cells("b(star label( )  fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS models. Dependent Variable: Change in Black Elected Officials (1964-1980)\label{baseline}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c c c} \toprule & (1) & (2) & (3) & (4) \\  \cmidrule(r){2-5}") ///
prefoot(" & & & & \\ State Trends & No  & Yes & Yes & Yes \\ Economic Controls & No & Yes & Yes & Yes \\ Other Controls & No & No & Yes & Yes \\ Coverage X Controls & No & No & No & Yes   \\ \midrule") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \setlength\labelsep{0pt} \item Robust standard errors in parenthesis clustered by judicial divisions. Controls: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}")



********************************************************************************
* TABLE 3
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_ShareBl_CountyGoverningBody ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2
xi: reg ch_ShareBl_JudLawEnf ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3
xi: reg ch_ShareBl_OtherCounty ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4
xi: reg ch_ShareBl_AllMunicipality ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5
xi: reg ch_ShareBl_AllEducation ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6

// Border sample
cd $dataset
use "dataset_border_1.dta", clear

twowayclusterBTF_withconst ch_ShareBl_CountyGoverningBody ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor2_clust_w
twowayclusterBTF_withconst ch_ShareBl_JudLawEnf ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor3_clust_w
twowayclusterBTF_withconst ch_ShareBl_OtherCounty ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor4_clust_w
twowayclusterBTF_withconst ch_ShareBl_AllMunicipality ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor5_clust_w
twowayclusterBTF_withconst ch_ShareBl_AllEducation ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor6_clust_w

cd $save_figures_tables
estout cme2 cme3 cme4 cme5 cme6 ///
using table_3.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n  0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) title(OLS models. Dependent Variable: Change in Black Elected Officials (1964-1980)\label{cme}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{adjustbox}{scale=0.85} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c c c c} \toprule & \multicolumn{3}{c}{County Governments} & \multicolumn{2}{c}{Other Governments} \\ \cmidrule(r){2-4} \cmidrule(r){5-6} & Commission & Judiciary and & Other & Municipality & School District\\ & &  Enforcement &  & & \\ & (1) & (2) & (3) & (4) & (5) \\ \midrule \multicolumn{6}{l}{\hspace{-0.3 cm} \emph{Panel A: Overall Sample}} \\ \midrule") ///
prefoot("& & & & &\\ Controls  & Yes & Yes & Yes & Yes & Yes\\ Controls X Coverage & Yes & Yes & Yes & Yes & Yes\\") ///
postfoot("\midrule \multicolumn{6}{l}{\hspace{-0.3 cm} \emph{Panel B: Border Sample}} \\ \midrule") ///

estout bor2_clust_w bor3_clust_w bor4_clust_w bor5_clust_w bor6_clust_w ///
using table_3.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
prefoot("& & & & & \\") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions in Panel A and by judicial divisions and border segments in Panel B. State trends are included in Panel A and county-pair and coverage trends are included in Panel B. Controls in Panel A: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes}\end{threeparttable} \end{adjustbox} \end{center} \end{table}")



********************************************************************************
* TABLE 4
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_AllCountyOfficials ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE ///
if ch_ShareBl_AllCountyOfficials!=. , ///
nocon robust cluster(judicial_divisions_id)
eststo den1
xi: reg ch_AllMunicipality ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE ///
if ch_ShareBl_AllMunicipality!=. , ///
nocon robust cluster(judicial_divisions_id)
eststo den2
xi: reg ch_AllEducation ///
black_share60_lit_nc black_share60 $interaction_full ///
ibn.STATE ///
if FIPSTATE!="51" & ch_ShareBl_AllEducation!=. & ch_ShareBl_AllEducation!=. , ///
nocon robust cluster(judicial_divisions_id)
eststo den3

cd $save_figures_tables
estout den1 den2 den3 ///
using table_4.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n 0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc black_share60) cells("b(star label( )  fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS models. Dependent Variable: Change in the number of members of elective bodies \label{den}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c c} \toprule & County & Municipality & School District\\ & (1) &(2) & (3) \\ \cmidrule(r){2-4}") ///
prefoot(" & & &  \\ Controls & Yes & Yes & Yes \\ Controls X Coverage & Yes & Yes & Yes \\ State Trends & Yes & Yes  & Yes \\ \midrule") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions. Controls: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}")



********************************************************************************
* TABLE 5
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

gen black_share60_AL = black_share60 * AL
gen black_share60_MIXED = black_share60 * MIXED

xi: reg ch_ShareBl_CountyGoverningBody ///
black_share60_lit_SMD black_share60_lit_nc ///
black_share60 black_share60_SMD ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo SMD82

xi: reg ch_ShareBl_CountyGoverningBody ///
black_share60_lit_SMD black_share60_lit_switchB ///
black_share60_lit_MMD_80 black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo SMD2

cd $save_figures_tables
estout SMD82 SMD2 using table_5.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_SMD black_share60_lit_nc black_share60_SMD black_share60 black_share60_lit_switchB black_share60_lit_MMD_80) order(black_share60_lit_SMD black_share60_lit_nc black_share60_SMD black_share60) cells("b(star label( )  fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS models. Dependent Variable: Change in Black Elected Officials, County Commission\label{SMD}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c} \toprule & (1) & (2) \\ \cmidrule(r){2-3}") ///
prefoot(" & & \\ Controls & Yes & Yes \\ Controls X Coverage & Yes & Yes \\ State Trends & Yes & Yes \\ \midrule") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions. Controls: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}")



********************************************************************************
* TABLE 6
********************************************************************************

cd $dataset
use "dataset_wide_4.dta", clear

xi: reg lndiff_county_current82 ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo spe2
xi: reg lndiff_county_capital82 ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo spe3
xi: reg lndiff_nocounty_current82 ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo nocounty2
xi: reg lndiff_nocounty_capital82 ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo nocounty3

// Border sample
cd $dataset
use "dataset_border_2.dta", clear

twowayclusterBTF_withconst lndiff_county_current_prepost ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
literacy_nc pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo spe2bor_clust_w
twowayclusterBTF_withconst lndiff_county_cap_prepost ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
literacy_nc pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo spe3bor_clust_w

cd $save_figures_tables
estout spe2 spe3 spe2bor_clust_w spe3bor_clust_w nocounty2 nocounty3 using table_6.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60_lit_SMD) order(black_share60_lit_SMD black_share60_lit_nc) cells("b(star label( )  fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS regressions. Dependent Variable: Change in Local Spending (1957-1982)\label{spending}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\footnotesize \begin{tabular}{@{}l c c c c c c} \toprule & \multicolumn{2}{c}{Overall Sample} & \multicolumn{2}{c}{Border Sample} & \multicolumn{2}{c}{Placebo}\\ \cmidrule(r){2-3} \cmidrule(r){4-5} \cmidrule(r){6-7} & (1) &(2) & (3) & (4) & (5) & (6) \\ & Current & Capital & Current & Capital & Current & Capital \\ \cmidrule(r){2-7}") ///
prefoot(" & & & & & & \\ Controls & Yes & Yes & Yes & Yes & Yes & Yes \\ Controls X Coverage & Yes & Yes & No & No & Yes & Yes \\ Controls X Coverage X SMD & Yes & Yes & No & No & Yes & Yes \\ \midrule") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions in columns (1)-(2)-(5)-(6) and by judicial divisions and border segments in columns (3)-(4). State trends are included in columns (1)-(2)-(5)-(6) and county-pair and coverage trends are included in columns (3)-(4). Controls in columns (1)-(2)-(5)-(6): \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}") 



********************************************************************************
* TABLE 7
********************************************************************************

cd $dataset
use "dataset_wide_4.dta", clear

xi: reg lndiff_county_totalIG82 ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo tax1
xi: reg lndiff_county_rev_own82 ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo tax2
xi: reg ch_share_fed_state ///
black_share60 black_share60_SMD black_share60_lit_nc black_share60_lit_SMD ///
$interaction_full_SMD i.FIPSTATE ///
if lndiff_county_totalIG82!=. , ///
nocon robust cluster(judicial_divisions_id) 
eststo share_fed_state

// Border sample
cd $dataset
use "dataset_border_2.dta", clear

twowayclusterBTF_withconst lndiff_county_totalIG82 ///
black_share60 black_share60_SMD black_share60_lit_nc ///
black_share60_lit_SMD literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo tax1bor_clust_w
twowayclusterBTF_withconst lndiff_county_rev_own82 ///
black_share60 black_share60_SMD black_share60_lit_nc ///
black_share60_lit_SMD literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo tax4bor_clust_w
twowayclusterBTF_withconst ch_share_fed_state ///
black_share60 black_share60_SMD black_share60_lit_nc ///
black_share60_lit_SMD literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 ///
if lndiff_county_totalIG82!=. , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo tax6bor_clust_w

cd $save_figures_tables
estout tax1 tax2 share_fed_state tax1bor_clust_w tax4bor_clust_w tax6bor_clust_w using table_7.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60_lit_SMD) order(black_share60_lit_SMD black_share60_lit_nc) cells("b(star label( )  fmt(4))" " se(par label( ) fmt(4))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS regressions. Dependent Variable: Change in Revenues (1957-1982)\label{revenues}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\footnotesize \begin{tabular}{@{}l c c c c c c} \toprule & \multicolumn{3}{c}{Overall Sample} & \multicolumn{3}{c}{Border Sample} \\ \cmidrule(r){2-4} \cmidrule(r){5-7} & (1) &(2) & (3) & (4) & (5) & (6) \\ & Intergovernmental & County Own & Federal & Intergovernmental & County Own  & Federal \\ & Revenues & Revenues & Share & Revenues & Revenues & Share \\ \cmidrule(r){2-7}") ///
prefoot(" & & & & & & \\ Controls & Yes & Yes & Yes & Yes & Yes & Yes \\ Controls X Coverage & Yes & Yes & Yes & No & No & No \\ Controls X Coverage X SMD & Yes & Yes & Yes & No & No & No \\ State Trends & Yes & Yes & Yes & No & No & No \\ County-Pair Trends & No & No & No & Yes & Yes & Yes \\ \midrule" ) ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions in columns (1)-(2)-(3) and by judicial divisions and border segments in columns (4)-(5)-(6). State trends are included in columns (1)-(2)-(3) and county-pair and coverage trends are included in columns (4)-(5)-(6). Controls in columns (1)-(2)-(3): \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}")


