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

Output files:   - Figures A.1 to A.5
				- Tables A.1 to A.7

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
global controls "urbanB60 unemp60 family_less_3000 pop60 school_low pro_black_county cotton_suitability cotton_share_land1964 anti_black_county rep_share_1964"
global beos "ShareBl_AllOfficials_1964 ShareBl_AllOfficials_1980 ShareBl_AllMunicipality_1964 ShareBl_AllMunicipality_1980 ShareBl_AllEducation_1964 ShareBl_AllEducation_1980 ShareBl_AllCountyOfficials_1964 ShareBl_AllCountyOfficials_1980 ShareBl_CountyGoverningBody_1964 ShareBl_CountyGoverningBody_1980 ShareBl_JudLawEnf_1964 ShareBl_JudLawEnf_1980 ShareBl_OtherCounty_1964 ShareBl_OtherCounty_1980"
global sumstats_controls "black_share60 population60 unemp60 urbanB60 family_less_3000 school_low cotton_suitability cotton_share_land1964 anti_black_county pro_black_county rep_share_1964"
global county_expenditure "county_expenditure57_re_pc county_expenditure82_re_pc county_current57_re_pc county_current82_re_pc county_cap57_re_pc county_cap82_re_pc"
global ln_county_expenditure "ln_county_expenditure57_real_pc ln_county_expenditure82_real_pc ln_county_current57_real_pc ln_county_current82_real_pc ln_county_cap57_real_pc ln_county_cap82_real_pc"
global interaction_full "c.urbanB60#literacy_nc c.unemp60#literacy_nc c.family_less_3000#literacy_nc c.pop60#literacy_nc c.school_low#literacy_nc c.cotton_suitability#literacy_nc c.cotton_share_land1964#literacy_nc c.anti_black_county#literacy_nc c.pro_black_county#literacy_nc c.rep_share_1964#literacy_nc"


********************************************************************************
* FIGURE A.1
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

spmap black_share60 using "countiescoord.dta", ///
id(geo) ndfcolor(gs14) legstyle(2) legjunction(" - ") legend(ring(1) ///
pos(6) size(*1.8)) title("") ndlabel("No data") ///
fcolor(white gs11 gs8 gs5 gs3 black) ///
clmethod(custom) clbreaks(0 10 20 30 40 50 100) legorder(lohi) ///
legend(label(2 "Below 10%") label(3 "10% to 20%") label(4 "20% to 30%") ///
label(5 "30% to 40%") label(6 "40% to 50%") label(7 "Above 50%") ///
color(dknavy*2.0))

cd $save_figures_tables
graph export "figure_a_1.eps", replace


********************************************************************************
* FIGURE A.2
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

destring FIPSTATE, replace
drop if literacy_nc == .

cd $build

// All officials
* Treatment
binsreg ch_ShareBl_AllOfficials ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("All Officials: Covered", size(medium))
graph save "binsreg_allofficials_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_AllOfficials ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("All Officials: Not Covered", size(medium)) 
graph save "binsreg_allofficials_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_allofficials_linearity_0_100_treatment.gph" ///
"binsreg_allofficials_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "all_linearity.gph", replace

// County governments
* Treatment
binsreg ch_ShareBl_CountyGoverningBody ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("County Gov. Bodies: Covered", size(medium)) 
graph save "binsreg_commissioners_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_CountyGoverningBody ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("County Gov. Bodies: Not Covered", size(medium)) 
graph save "binsreg_commissioners_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_commissioners_linearity_0_100_treatment.gph" ///
"binsreg_commissioners_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "commission_linearity.gph", replace

// School district
* Treatment
binsreg ch_ShareBl_AllEducation ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("School District: Covered", size(medium)) 
graph save "binsreg_education_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_AllEducation ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("School District: Not Covered", size(medium)) 
graph save "binsreg_education_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_education_linearity_0_100_treatment.gph" ///
"binsreg_education_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "education_linearity.gph", replace

// Municipality
* Treatment
binsreg ch_ShareBl_AllMunicipality ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Municipality: Covered", size(medium)) 
graph save "binsreg_municipality_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_AllMunicipality ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Municipality: Not Covered", size(medium)) 
graph save "binsreg_municipality_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_municipality_linearity_0_100_treatment.gph" ///
"binsreg_municipality_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "municipality_linearity.gph", replace

// Judiciary and law enforcement
* Treatment
binsreg ch_ShareBl_JudLawEnf ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Judiciary and Enforc.: Covered", size(medium)) 
graph save "binsreg_judlaw_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_JudLawEnf ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Judiciary and Enforc.: Not Covered", size(medium)) 
graph save "binsreg_judlaw_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_judlaw_linearity_0_100_treatment.gph" ///
"binsreg_judlaw_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "judlaw_linearity.gph", replace

// Other officials
* Treatment
binsreg ch_ShareBl_OtherCounty ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Other: Covered", size(medium)) 
graph save "binsreg_othercounty_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_OtherCounty ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Other: Not Covered", size(medium)) 
graph save "binsreg_othercounty_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_othercounty_linearity_0_100_treatment.gph" ///
"binsreg_othercounty_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "othercounty_linearity.gph", replace

graph combine "all_linearity.gph" "commission_linearity.gph" ///
"judlaw_linearity.gph" "othercounty_linearity.gph" ///
"education_linearity.gph" "municipality_linearity.gph", ///
graphregion(color(white)) rows(3) ycommon

cd $save_figures_tables
graph export "figure_a_2.eps", replace



********************************************************************************
* FIGURE A.3
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

destring FIPSTATE, replace
drop if literacy_nc == .

cd $build

// Common support
keep if black_share60 < 68.9

// All officials
* Treatment
binsreg ch_ShareBl_AllOfficials ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("All Officials: Covered", size(medium))
graph save "binsreg_allofficials_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_AllOfficials ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("All Officials: Not Covered", size(medium)) 
graph save "binsreg_allofficials_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_allofficials_linearity_0_100_treatment.gph" ///
"binsreg_allofficials_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "all_linearity_cs.gph", replace

// County governments
* Treatment
binsreg ch_ShareBl_CountyGoverningBody ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("County Gov. Bodies: Covered", size(medium)) 
graph save "binsreg_commissioners_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_CountyGoverningBody ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("County Gov. Bodies: Not Covered", size(medium)) 
graph save "binsreg_commissioners_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_commissioners_linearity_0_100_treatment.gph" ///
"binsreg_commissioners_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "commission_linearity_cs.gph", replace

// School district
* Treatment
binsreg ch_ShareBl_AllEducation ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("School District: Covered", size(medium)) 
graph save "binsreg_education_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_AllEducation ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("School District: Not Covered", size(medium)) 
graph save "binsreg_education_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_education_linearity_0_100_treatment.gph" ///
"binsreg_education_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "education_linearity_cs.gph", replace

// Municipality
* Treatment
binsreg ch_ShareBl_AllMunicipality ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Municipality: Covered", size(medium)) 
graph save "binsreg_municipality_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_AllMunicipality ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Municipality: Not Covered", size(medium)) 
graph save "binsreg_municipality_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_municipality_linearity_0_100_treatment.gph" ///
"binsreg_municipality_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "municipality_linearity_cs.gph", replace

// Judiciary and law enforcement
* Treatment
binsreg ch_ShareBl_JudLawEnf ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Judiciary and Enforc.: Covered", size(medium)) 
graph save "binsreg_judlaw_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_JudLawEnf ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Judiciary and Enforc.: Not Covered", size(medium)) 
graph save "binsreg_judlaw_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_judlaw_linearity_0_100_treatment.gph" ///
"binsreg_judlaw_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "judlaw_linearity_cs.gph", replace

// Other officials
* Treatment
binsreg ch_ShareBl_OtherCounty ///
black_share60 $controls i.FIPSTATE if literacy_nc == 1, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(O O) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Other: Covered", size(medium)) 
graph save "binsreg_othercounty_linearity_0_100_treatment.gph", replace
* Control
binsreg ch_ShareBl_OtherCounty ///
black_share60 $controls i.FIPSTATE if literacy_nc == 0, ///
by(literacy_nc) vce(cluster judicial_divisions_id) ///
nbins(5) samebinsby binspos(es) dotsgrid(1) cb(1,1) ///
cbplotopt(fcolor(gs13) lcolor(white) lwidth(thin)) polyreg() ///
graphregion(color(white)) bycolors(black) bysymbols(T T) legend(off) ///
ylabel( , glwidth(0.00) glpattern(shortdash) glcolor(erose)) ///
xscale(range(0 100)) xlabel(0(20)100) ytitle("") xtitle("") ///
title("Other: Not Covered", size(medium)) 
graph save "binsreg_othercounty_linearity_0_100_control.gph", replace
* Combined
graph combine "binsreg_othercounty_linearity_0_100_treatment.gph" ///
"binsreg_othercounty_linearity_0_100_control.gph", ///
graphregion(color(white)) ycommon
graph save "othercounty_linearity_cs.gph", replace

// Combine all the plots
graph combine "all_linearity_cs.gph" "commission_linearity_cs.gph" ///
"judlaw_linearity_cs.gph" "othercounty_linearity_cs.gph" ///
"education_linearity_cs.gph" "municipality_linearity_cs.gph",  ///
graphregion(color(white)) rows(3)

cd $save_figures_tables
graph export "figure_a_3.eps", replace



********************************************************************************
* FIGURE A.4
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

gen literacy_nc_SMD = literacy_nc * SMD

foreach var of varlist black_share60 black_share40 black_share50 pop60 pop50 ///
unemp60 unemp50 family_less_3000 family_less50_2000 urbanB60 urban50 ///
school_low school_low50 cotton_suitability pro_black_county anti_black_county {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc literacy_nc_SMD SMD , ///
robust cluster(judicial_divisions_id)
est sto `var'_std

}

coefplot (black_share60_std, label(Black Share) msymbol(O) msize(medium) ///
color(black)ciopts(lcolor(black))) ///
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
, keep(literacy_nc_SMD) title(SMD: Overall Sample) levels(95) ytitle("") ///
ylabel("") legend(rows(3) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))
cd $build
graph save "figure_a_4_a.gph", replace

// Border sample
cd $dataset
use "dataset_border_1.dta", clear

foreach var of varlist black_share60 black_share40 black_share50 pop60 pop50 ///
unemp60 unemp50 family_less_3000 family_less50_2000 urbanB60 urban50 ///
school_low school_low50 cotton_suitability pro_black_county anti_black_county {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc literacy_nc_SMD SMD , ///
robust cluster(judicial_border)
est sto `var'_std

}

coefplot (black_share60_std, label(Black Share) msymbol(O) msize(medium) ///
color(black)ciopts(lcolor(black))) ///
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
, keep(literacy_nc_SMD) title(SMD: Border Sample) levels(95) ytitle("") ///
ylabel("") legend(rows(3) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))
cd $build
graph save "figure_a_4_b.gph", replace

// Combine all the plots
cd $build
graph combine "figure_a_4_a.gph" "figure_a_4_b.gph", ///
graphregion(color(white)) xcommon ycommon

cd $save_figures_tables
graph export "figure_a_4.eps", replace



********************************************************************************
* FIGURE A.5
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

gen literacy_nc_SMD = literacy_nc * SMD

foreach var of varlist unempchange60_50 ruralchange60_50 urbanchange60_50 ///
familypovchange60_50 diffshareblack60_50 lndiff_pop_50_60 ///
school_lowchange60_50 {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc literacy_nc_SMD SMD , ///
robust cluster(judicial_divisions_id)
est sto `var'_std 

}

coefplot (diffshareblack60_50_std, label({&Delta} Black share) msymbol(O) ///
msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange60_50_std, label({&Delta} Unemployment) msymbol(Oh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange60_50_std, label({&Delta} Poverty) msymbol(S) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange60_50_std, label({&Delta} Unskilled) msymbol(Sh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange60_50_std, label({&Delta} Percent Urban) msymbol(D) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_50_60_std , label({&Delta} Population)msymbol(Dh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc_SMD) title(Pre-VRA SMD: Overall) levels (95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))
cd $build
graph save "figure_a_5_a.gph", replace

foreach var of varlist unempchange ruralchange urbanchange familypovchange ///
diffshareblack80 lndiff_pop_60_80 school_lowchange {

quietly summarize `var'
gen `var'_std = (`var'-r(mean)) / r(sd)
xi: reg `var'_std literacy_nc literacy_nc_SMD SMD, ///
robust cluster(judicial_divisions_id)
est sto `var'_std

}

coefplot (diffshareblack80_std, label({&Delta} Black share) msymbol(O) ///
msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange_std, label({&Delta} Unemployment) msymbol(Oh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange_std, label({&Delta} Poverty) msymbol(S) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange_std, label({&Delta} Unskilled) msymbol(Sh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange_std, label({&Delta} Percent Urban) msymbol(D) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_60_80_std , label({&Delta} Population) msymbol(Dh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc_SMD) title(Post-VRA SMD: Overall) levels (95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))
cd $build
graph save "figure_a_5_b.gph", replace


// Border sample
cd $dataset
use "dataset_border_1.dta", clear

foreach var of varlist unempchange60_50 ruralchange60_50 urbanchange60_50 ///
familypovchange60_50 diffshareblack60_50 lndiff_pop_50_60 ///
school_lowchange60_50 {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc literacy_nc_SMD SMD , ///
robust cluster(judicial_border)
est sto `var'_std

}

coefplot (diffshareblack60_50_std, label({&Delta} Black share) msymbol(O) ///
msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange60_50_std, label({&Delta} Unemployment) msymbol(Oh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange60_50_std, label({&Delta} Poverty) msymbol(S) ///
msize(medium) color(black) ciopts(lcolor(black)))  ///
(school_lowchange60_50_std, label({&Delta} Unskilled) msymbol(Sh) ///
msize(medium) color(black) ciopts(lcolor(black)))    ///
(urbanchange60_50_std, label({&Delta} Percent Urban) msymbol(D) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_50_60_std, label({&Delta} Population)msymbol(Dh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc_SMD) title(Pre-VRA SMD: Border)  levels (95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))
cd $build
graph save "figure_a_5_c.gph", replace


foreach var of varlist unempchange ruralchange urbanchange familypovchange ///
diffshareblack80 lndiff_pop_60_80 school_lowchange {

quietly summarize `var'
gen `var'_std = (`var' - r(mean)) / r(sd)
xi: reg `var'_std literacy_nc literacy_nc_SMD SMD , ///
robust cluster(judicial_border)
est sto `var'_std

}

coefplot (diffshareblack80_std, label({&Delta} Black share) msymbol(O) ///
msize(medium) color(black)ciopts(lcolor(black))) ///
(unempchange_std, label({&Delta} Unemployment) msymbol(Oh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(familypovchange_std, label({&Delta} Poverty) msymbol(S) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(school_lowchange_std, label({&Delta} Unskilled) msymbol(Sh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(urbanchange_std, label({&Delta} Percent Urban) msymbol(D) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
(lndiff_pop_60_80_std, label({&Delta} Population) msymbol(Dh) ///
msize(medium) color(black) ciopts(lcolor(black))) ///
, keep(literacy_nc_SMD) title(Post-VRA SMD: Border) levels (95) ytitle("") ///
ylabel("") legend(rows(2) size(small)) xline(0, lpattern(dash)) ///
xlabels( , labs(medium) labgap() format(%9.1f))
cd $build
graph save "figure_a_5_d.gph", replace

// Combine all the plots
cd $build
graph combine "figure_a_5_a.gph" "figure_a_5_b.gph" "figure_a_5_c.gph" "figure_a_5_d.gph", ///
graphregion(color(white)) xcommon ycommon

cd $save_figures_tables
graph export "figure_a_5.eps", replace



********************************************************************************
* TABLE A.1
********************************************************************************

cd $dataset
use "dataset_wide_5.dta", clear

xi: estpost summarize $beos if literacy_nc == 1
eststo sum1
xi: estpost summarize $beos if literacy_nc == 0
eststo sum2

xi: estpost summarize $sumstats_controls if literacy_nc == 1 
eststo sum9
xi: estpost summarize $sumstats_controls if literacy_nc == 0
eststo sum10

xi: estpost summarize $county_expenditure if literacy_nc == 1 
eststo sum5a
xi: estpost summarize $county_expenditure if literacy_nc == 0 
eststo sum6a

xi: estpost summarize $ln_county_expenditure if literacy_nc == 1 
eststo sum5b
xi: estpost summarize $ln_county_expenditure if literacy_nc == 0 
eststo sum6b

cd $save_figures_tables
estout sum1 sum2 using table_a_1.tex, replace cells("mean(fmt(2)) sd(fmt(2))") mlabels(none) collabels(none) label varwidth(12) modelwidth(10) style(tex) ///
title(Summary Statistics\label{summary_stat_new_jpe}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\scriptsize \begin{tabular}{@{}l c c c c} \toprule &  \multicolumn{2}{c}{Covered} & \multicolumn{2}{c}{Not covered} \\ \cmidrule(r){2-3} \cmidrule(r){4-5} & Mean & St. Dev. & Mean & St. Dev. \\ \midrule \rowgroup{\emph{Local Black Elected Officials}} \\ & & & & \\")

estout sum9 sum10 using table_a_1.tex, cells("mean(fmt(2)) sd(fmt(2))") append mlabels(none) collabels(none) label varwidth(12) modelwidth(10) style(tex) ///
posthead("\midrule \rowgroup{\emph{County Characteristics}} \\ & & & & \\") ///
postfoot(" & & & & \\ Counties & \multicolumn{2}{c}{511} & \multicolumn{2}{c}{537} \\ \bottomrule \end{tabular} \begin{tablenotes} \item \end{tablenotes} \end{threeparttable} \end{center} \end{table}")



********************************************************************************
* TABLE A.2
********************************************************************************

cd $save_figures_tables
estout sum5a sum6a using table_a_2.tex, replace cells("mean(fmt(2)) sd(fmt(2))") mlabels(none) collabels(none) label varwidth(12) modelwidth(10) style(tex) ///
title(Summary Statistics\label{summary_stat2final_new_jpe}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\scriptsize \begin{tabular}{@{}l c c c c} \toprule & \multicolumn{2}{c}{Covered} & \multicolumn{2}{c}{Not covered} \\ \cmidrule(r){2-3} \cmidrule(r){4-5} & Mean & St. Dev. & Mean & St. Dev. \\ \midrule \rowgroup{\emph{Local Expenditure, percapita (real 2000 dollars)}} \\ & & & & \\ ")

estout sum5b sum6b using table_a_2.tex, cells("mean(fmt(2)) sd(fmt(2))") append mlabels(none) collabels(none) label varwidth(12) modelwidth(10) style(tex) ///
posthead("\midrule \rowgroup{\emph{ln Local Expenditure, percapita (real 2000 dollars)}} \\ & & & & \\") ///
postfoot(" & & & & \\ \bottomrule \end{tabular} \begin{tablenotes} \item \end{tablenotes} \end{threeparttable} \end{center} \end{table}")



********************************************************************************
* TABLE A.3
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if black_share60<68.9 , ///
nocon robust cluster(judicial_divisions_id)
eststo cs1
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if black_share60<68.9 , ///
nocon robust cluster(judicial_divisions_id)
eststo cs3
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if black_share60<68.9 , ///
nocon robust cluster(judicial_divisions_id)
eststo cs4
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if black_share60<68.9 , ///
nocon robust cluster(judicial_divisions_id)
eststo cs5
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if black_share60<68.9 , ///
nocon robust cluster(judicial_divisions_id)
eststo cs6

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE [aw=pop60] , ///
nocon robust cluster(judicial_divisions_id)
eststo pop1
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE [aw=pop60] , ///
nocon robust cluster(judicial_divisions_id)
eststo pop3
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE [aw=pop60] , ///
nocon robust cluster(judicial_divisions_id)
eststo pop4
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE [aw=pop60] , ///
nocon robust cluster(judicial_divisions_id)
eststo pop5
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE [aw=pop60] , ///
nocon robust cluster(judicial_divisions_id)
eststo pop6

cd $save_figures_tables
estout cs1 cs3 cs4 cs5 cs6 using table_a_3.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n  0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS models. Common Support and Population weighting. Dependent Variable: Change in Black Elected Officials (1964-1980)\label{cs_pop}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{adjustbox}{scale=0.85} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}lc c c c c} \toprule & \multicolumn{3}{c}{County Governments} & \multicolumn{2}{c}{Other Governments}\\ \cmidrule(r){2-4} \cmidrule(r){5-6} & Commission &  Judiciary and & Other & Municipality & School District\\ & &  Enforcement & & & \\ & (1) & (2) & (3) & (4) & (5) \\ \midrule \multicolumn{5}{l}{\hspace{-0.3 cm} \emph{Panel A: Common Support}} \\ \midrule") ///
prefoot("& & & & & \\ Controls & Yes & Yes & Yes & Yes & Yes \\ Controls X Coverage & Yes & Yes & Yes & Yes & Yes \\") ///
postfoot("\midrule \multicolumn{5}{l}{\hspace{-0.3 cm} \emph{Panel B: Population weighted}} \\ \midrule") ///

estout pop1 pop3 pop4 pop5 pop6 using table_a_3.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
prefoot("& & & & & \\ Controls & Yes & Yes & Yes & Yes & Yes \\ Controls X Coverage & Yes & Yes & Yes & Yes & Yes \\") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions. State trends are included. Controls: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes}\end{threeparttable} \end{adjustbox} \end{center} \end{table}")



********************************************************************************
* TABLE A.4
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$controls ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_noint
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$controls ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_noint
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$controls ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_noint
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$controls ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_noint
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$controls ibn.STATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_noint

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$interaction_full ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$interaction_full ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$interaction_full ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$interaction_full ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$interaction_full ibn.STATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6

// End at 1976
drop ch_ShareBl_AllOfficials
gen ch_ShareBl_AllOfficials = ShareBl_AllOfficials_1976 - ShareBl_AllOfficials_1964
drop ch_ShareBl_CountyGoverningBody
gen ch_ShareBl_CountyGoverningBody = ShareBl_CountyGoverningBody_1976 - ShareBl_CountyGoverningBody_1964
drop ch_ShareBl_Mun_Educ_Total
gen ch_ShareBl_Mun_Educ_Total = ShareBl_Mun_Educ_Total_1976 - ShareBl_Mun_Educ_Total_1964
drop ch_ShareBl_CountyGoverningBody
gen ch_ShareBl_CountyGoverningBody = ShareBl_CountyGoverningBody_1976 - ShareBl_CountyGoverningBody_1964
drop ch_ShareBl_JudLawEnf
gen ch_ShareBl_JudLawEnf = ShareBl_JudLawEnf_1976 - ShareBl_JudLawEnf_1964
drop ch_ShareBl_OtherCounty
gen ch_ShareBl_OtherCounty = ShareBl_OtherCounty_1976 - ShareBl_OtherCounty_1964
drop ch_ShareBl_AllMunicipality
gen ch_ShareBl_AllMunicipality = ShareBl_AllMunicipality_1976 - ShareBl_AllMunicipality_1964
drop ch_ShareBl_AllEducation
gen ch_ShareBl_AllEducation = ShareBl_AllEducation_1976 - ShareBl_AllEducation_1964

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_76_noint
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_76_noint
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_76_noint
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_76_noint
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_76_noint	

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_76
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_76
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_76
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_76
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_76		

// End at 1978
drop ch_ShareBl_AllOfficials
gen ch_ShareBl_AllOfficials = ShareBl_AllOfficials_1978 - ShareBl_AllOfficials_1964
drop ch_ShareBl_CountyGoverningBody
gen ch_ShareBl_CountyGoverningBody = ShareBl_CountyGoverningBody_1978 - ShareBl_CountyGoverningBody_1964
drop ch_ShareBl_Mun_Educ_Total
gen ch_ShareBl_Mun_Educ_Total = ShareBl_Mun_Educ_Total_1978 - ShareBl_Mun_Educ_Total_1964
drop ch_ShareBl_CountyGoverningBody
gen ch_ShareBl_CountyGoverningBody = ShareBl_CountyGoverningBody_1978 - ShareBl_CountyGoverningBody_1964
drop ch_ShareBl_JudLawEnf
gen ch_ShareBl_JudLawEnf = ShareBl_JudLawEnf_1978 - ShareBl_JudLawEnf_1964
drop ch_ShareBl_OtherCounty
gen ch_ShareBl_OtherCounty = ShareBl_OtherCounty_1978 - ShareBl_OtherCounty_1964
drop ch_ShareBl_AllMunicipality
gen ch_ShareBl_AllMunicipality = ShareBl_AllMunicipality_1978 - ShareBl_AllMunicipality_1964
drop ch_ShareBl_AllEducation
gen ch_ShareBl_AllEducation = ShareBl_AllEducation_1978 - ShareBl_AllEducation_1964	

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_78_noint
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_78_noint
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_78_noint
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_78_noint
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_78_noint		

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_78
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_78
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_78
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_78
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_78

// End at 1978-80
egen aver_ShareBl_AllOfficials = rowmean(ShareBl_AllOfficials_1978 ShareBl_AllOfficials_1980)
drop ch_ShareBl_AllOfficials
gen ch_ShareBl_AllOfficials = aver_ShareBl_AllOfficials - ShareBl_AllOfficials_1964
egen aver_ShareBl_CountyGoverningBody = rowmean(ShareBl_CountyGoverningBody_1978 ShareBl_CountyGoverningBody_1980)
egen aver_ShareBl_Mun_Educ_Total = rowmean(ShareBl_Mun_Educ_Total_1978 ShareBl_Mun_Educ_Total_1980)
drop ch_ShareBl_CountyGoverningBody
gen ch_ShareBl_CountyGoverningBody = aver_ShareBl_CountyGoverningBody - ShareBl_CountyGoverningBody_1964
drop ch_ShareBl_Mun_Educ_Total
gen ch_ShareBl_Mun_Educ_Total = aver_ShareBl_Mun_Educ_Total - ShareBl_Mun_Educ_Total_1964
egen aver_ShareBl_JudLawEnf = rowmean(ShareBl_JudLawEnf_1978 ShareBl_JudLawEnf_1980)
egen aver_ShareBl_OtherCounty = rowmean(ShareBl_OtherCounty_1978 ShareBl_OtherCounty_1980)
egen aver_ShareBl_AllMunicipality = rowmean(ShareBl_AllMunicipality_1978 ShareBl_AllMunicipality_1980)
egen aver_ShareBl_AllEducation = rowmean(ShareBl_AllEducation_1978 ShareBl_AllEducation_1980)
drop ch_ShareBl_JudLawEnf
gen ch_ShareBl_JudLawEnf = aver_ShareBl_JudLawEnf - ShareBl_JudLawEnf_1964
drop ch_ShareBl_OtherCounty
gen ch_ShareBl_OtherCounty = aver_ShareBl_OtherCounty - ShareBl_OtherCounty_1964
drop ch_ShareBl_AllMunicipality
gen ch_ShareBl_AllMunicipality = aver_ShareBl_AllMunicipality - ShareBl_AllMunicipality_1964
drop ch_ShareBl_AllEducation
gen ch_ShareBl_AllEducation = aver_ShareBl_AllEducation - ShareBl_AllEducation_1964

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme2_78_80_noint
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme3_78_80_noint
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme4_78_80_noint
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme5_78_80_noint
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions)
eststo cme6_78_80_noint				

xi: reg ch_ShareBl_CountyGoverningBody black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme2_78_80
xi: reg ch_ShareBl_JudLawEnf black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme3_78_80
xi: reg ch_ShareBl_OtherCounty black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme4_78_80
xi: reg ch_ShareBl_AllMunicipality black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions)
eststo cme5_78_80
xi: reg ch_ShareBl_AllEducation black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions)
eststo cme6_78_80	

xi: reg ch_ShareBl_CountyGoverningBody90 black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_90
xi: reg ch_ShareBl_JudLawEnf90 black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_90
xi: reg ch_ShareBl_OtherCounty90 black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_90
xi: reg ch_ShareBl_AllMunicipality90 black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_90
xi: reg ch_ShareBl_AllEducation90 black_share60_lit_nc black_share60 ///
$interaction_full i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_90

xi: reg ch_ShareBl_CountyGoverningBody90 black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme2_90_noint
xi: reg ch_ShareBl_JudLawEnf90 black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme3_90_noint
xi: reg ch_ShareBl_OtherCounty90 black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme4_90_noint
xi: reg ch_ShareBl_AllMunicipality90 black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE , ///
nocon robust cluster(judicial_divisions_id)
eststo cme5_90_noint
xi: reg ch_ShareBl_AllEducation90 black_share60_lit_nc black_share60 ///
$controls i.FIPSTATE if FIPSTATE!="51" , ///
nocon robust cluster(judicial_divisions_id)
eststo cme6_90_noint

cd $save_figures_tables
estout cme2_noint cme3_noint cme4_noint cme5_noint cme6_noint cme2 cme3 cme4 cme5 cme6 using table_a_4.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n  0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc) order(black_share60_lit_nc) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
title(OLS models. Dependent Variable: Change in Black Elected Officials for Different Periods\label{BEO_1976}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{adjustbox}{scale=0.7} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c c c c c c c c c} \toprule & \multicolumn{5}{c}{\emph{Without Coverage X Controls} } & \multicolumn{5}{c}{\emph{With Coverage X Controls}} \\ \cmidrule(r){2-6} \cmidrule(r){7-11} & Commission &  Judiciary & Other & Municipality & School District & Commission & Judiciary & Other & Municipality & School District \\ \\ & (1) & (2) & (3) & (4) & (5) & (6) & (7) & (8) & (9) & (10) \\ \midrule \multicolumn{6}{l}{\hspace{-0.3 cm} \emph{Panel A: Benchmark (1964-1980)}} \\ \midrule") ///
prefoot("& & & & & & & & & &\\") ///
postfoot( "\midrule \multicolumn{6}{l}{\hspace{-0.3 cm} \emph{Panel B: 1964-1976 }} \\ \midrule") ///
			
estout cme2_76_noint cme3_76_noint cme4_76_noint cme5_76_noint cme6_76_noint cme2_76 cme3_76 cme4_76 cme5_76 cme6_76 using table_a_4.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc) order(black_share60_lit_nc) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
prefoot("& & & & & & & & & &\\") ///
postfoot( "\midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel C: 1964-1978 }} \\ \midrule") ///
			
estout cme2_78_noint cme3_78_noint cme4_78_noint cme5_78_noint cme6_78_noint cme2_78 cme3_78 cme4_78 cme5_78 cme6_78 using table_a_4.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc) order(black_share60_lit_nc) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
prefoot("& & & & & & & & & &\\") ///
postfoot( "\midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel D: 1964- Average 1978-80 }} \\ \midrule") ///
						
estout cme2_78_80_noint cme3_78_80_noint cme4_78_80_noint cme5_78_80_noint cme6_78_80_noint cme2_78_80 cme3_78_80 cme4_78_80 cme5_78_80 cme6_78_80 using table_a_4.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc) order(black_share60_lit_nc) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
prefoot("& & & & & & & & & &\\") ///
postfoot( "\midrule \multicolumn{4}{l}{\hspace{-0.3 cm} \emph{Panel E: 1964-1990 }} \\ \midrule") ///

estout cme2_90_noint cme3_90_noint cme4_90_noint cme5_90_noint cme6_90_noint cme2_90 cme3_90 cme4_90 cme5_90 cme6_90 using table_a_4.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc) order(black_share60_lit_nc) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
prefoot("& & & & & & & & & &\\") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions. State trends are included. Controls: \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes}\end{threeparttable} \end{adjustbox} \end{center} \end{table}")



********************************************************************************
* TABLE A.5
********************************************************************************

cd $dataset
use "dataset_border_1.dta", clear

twowayclusterBTF_withconst ch_ShareBl_CountyGoverningBody ///
black_share60 black_share60_lit_nc rep_share_1964 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor2_clust_w_2
twowayclusterBTF_withconst ch_ShareBl_JudLawEnf ///
black_share60 black_share60_lit_nc rep_share_1964 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor3_clust_w_2
twowayclusterBTF_withconst ch_ShareBl_OtherCounty ///
black_share60 black_share60_lit_nc rep_share_1964 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor4_clust_w_2
twowayclusterBTF_withconst ch_ShareBl_AllMunicipality ///
black_share60 black_share60_lit_nc rep_share_1964 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor5_clust_w_2
twowayclusterBTF_withconst ch_ShareBl_AllEducation ///
black_share60 black_share60_lit_nc rep_share_1964 literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor6_clust_w_2

twowayclusternoweightsBTF ch_ShareBl_CountyGoverningBody ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor2_clust
twowayclusternoweightsBTF ch_ShareBl_JudLawEnf ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor3_clust
twowayclusternoweightsBTF ch_ShareBl_OtherCounty ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor4_clust
twowayclusternoweightsBTF ch_ShareBl_AllMunicipality ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor5_clust
twowayclusternoweightsBTF ch_ShareBl_AllEducation ///
black_share60 black_share60_lit_nc literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor6_clust

cd $save_figures_tables
estout bor2_clust_w_2 bor3_clust_w_2 bor4_clust_w_2 bor5_clust_w_2 bor6_clust_w_2 using table_a_5.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) rename (0.literacy_nc#c.unemp60 unemp60_n 0.literacy_nc#c.pop60 pop60_n 0.literacy_nc#c.family_less_3000 family_less_3000_n 0.literacy_nc#c.school_low school_low_n 0.literacy_nc#c.urbanB60 urbanB60_n  0.literacy_nc#c.pro_black_county pro_black_county_n 0.literacy_nc#c.cotton_suitability cotton_suitability_n 0.literacy_nc#c.anti_black_county anti_black_county_n) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
title(OLS models. Dependent Variable: Change in Black Elected Officials (1964-1980)\label{bor_rob}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{adjustbox}{scale=0.85} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1.2}\small \begin{tabular}{@{}l c c c c c} \toprule & \multicolumn{3}{c}{County Governments} & \multicolumn{2}{c}{Other Governments} \\ \cmidrule(r){2-4} \cmidrule(r){5-6} & Commission &  Judiciary and & Other & Municipality & School District\\ & & Enforcement &  & & \\ & (1) & (2) & (3) & (4) & (5) \\ \midrule \multicolumn{6}{l}{\hspace{-0.3 cm} \emph{Panel A: Republican share control}} \\ \midrule")  ///
prefoot("& & & & &\\") ///
postfoot( "\midrule \multicolumn{6}{l}{\hspace{-0.3 cm} \emph{Panel B: Without border weights}} \\ \midrule") ///

estout bor2_clust bor3_clust bor4_clust bor5_clust bor6_clust using table_a_5.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60) order(black_share60_lit_nc black_share60) cells("b(star label( ) fmt(3))" " se(par label( ) fmt(3))") append mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
prefoot("& & & & &\\") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions and border segments. County-pair and coverage trends are included. ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes}\end{threeparttable} \end{adjustbox} \end{center} \end{table}")



********************************************************************************
* TABLE A.6
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_ShareBl_CountyGoverningBody ///
black_share60_lit_SMD black_share60_lit_nc black_share60 black_share60_SMD ///
$interaction_full i.FIPSTATE [aw=pop60] , ///
nocon robust cluster(judicial_divisions_id)
eststo pop2_SMD
xi: reg ch_ShareBl_CountyGoverningBody ///
black_share60_lit_SMD black_share60_lit_nc black_share60 black_share60_SMD ///
$interaction_full i.FIPSTATE if black_share60<68.9 , ///
nocon robust  cluster(judicial_divisions_id)
eststo cs2_SMD

cd $dataset
use "dataset_border_1.dta", clear

twowayclusterBTF_withconst ch_ShareBl_CountyGoverningBody ///
black_share60_lit_SMD black_share60_lit_nc ///
black_share60 black_share60_SMD literacy_nc ///
pair_fixedeffect_1-pair_fixedeffect_127 , ///
w(county_w) fcluster(judicial_border) tcluster(judicial_divisions_id)
eststo bor8_clust_w

cd $dataset
use "dataset_wide_1.dta", clear

gen black_share60_AL = black_share60 * AL
gen black_share60_MIXED = black_share60 * MIXED

xi: reg ch_ShareBl_Mun_Educ_Total black_share60_lit_SMD black_share60_lit_nc black_share60 black_share60_SMD $interaction_full i.FIPSTATE , nocon robust cluster(judicial_divisions_id)
eststo SMD_munedu

cd $save_figures_tables
estout pop2_SMD cs2_SMD bor8_clust_w SMD_munedu using table_a_6.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_SMD black_share60_lit_nc) order(black_share60_lit_SMD black_share60_lit_nc) cells("b(star label( )  fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth (10) style(tex) ///
title(OLS models. Dependent Variable: Change in Black Elected Officials, by County Commission Election Rule \label{SMD_robust}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{1}\small \begin{tabular}{@{}l c c c c} \toprule & \multicolumn{3}{c}{County Commission} & \multicolumn{1}{c}{Municipality and School District} \\ \cmidrule(r){2-4} \cmidrule(r){5-5} & (1) &(2) & (3) & (4) \\ & Population Weighted & Common Support &  Border Sample \\ \cmidrule(r){2-4} \cmidrule(r){5-5}") ///
prefoot(" & & & & \\ Controls & Yes & Yes & Yes & Yes \\ Controls X Coverage & Yes  & Yes & No & Yes \\ \midrule" ) ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \item Robust standard errors in parenthesis clustered by judicial divisions in columns (1)-(2)-(4) and by judicial divisions and border segments in column (3). State trends are included in columns (1)-(2)-(4) and county-pair and coverage trends are included in column (3). Controls in columns (1)-(2)-(4): \emph{Unemployment rate (\%), 1960; Families below poverty line (\%), 1960; Low-skilled (\%), 1960; Population, 1960; Urban population (\%), 1960; Agricultural productivity; Cotton share (\%), 1964; Pro-Black protest, 1960-64; Anti-Black protest, 1960-64; Republican share (\%), 1964.} ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}")



********************************************************************************
* TABLE A.7
********************************************************************************

cd $dataset
use "dataset_wide_1.dta", clear

xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 , ///
nocon robust cluster(judicial_divisions_id)
eststo base1_c
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 unemp60 family_less_3000 pop60 ///
school_low urbanB60 cotton_suitability cotton_share_land1964 ibn.STATE , ///
nocon robust  cluster(judicial_divisions_id)
eststo base2_c
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 unemp60 pop60 family_less_3000 ///
school_low urbanB60 cotton_suitability cotton_share_land1964 ///
pro_black_county anti_black_county rep_share_1964 ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo base3_c
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 urbanB60 unemp60 family_less_3000 ///
pop60 school_low pro_black_county cotton_suitability cotton_share_land1964 ///
anti_black_county rep_share_1964 unemp60_lit family_less_3000_lit ///
pop60_lit school_low_lit urbanB60_lit  cotton_suitability_lit ///
cotton_share_land1964_lit pro_black_county_lit anti_black_county_lit ///
rep_share_1964_lit ibn.STATE , ///
nocon robust cluster(judicial_divisions_id)
eststo base4_c
xi: reg ch_ShareBl_AllOfficials ///
black_share60_lit_nc black_share60 urbanB60 unemp60 family_less_3000 ///
pop60 school_low pro_black_county cotton_suitability cotton_share_land1964 ///
anti_black_county rep_share_1964 unemp60_lit family_less_3000_lit ///
pop60_lit school_low_lit urbanB60_lit  cotton_suitability_lit ///
cotton_share_land1964_lit pro_black_county_lit anti_black_county_lit ///
rep_share_1964_lit ibn.STATE , ///
nocon robust cluster(county)
eststo base5_c

cd $save_figures_tables
estout base1_c base2_c base3_c base4_c base5_c using table_a_7.tex, stats(r2_a  N, fmt(3 0) labels("Adj. R-Square" "N")) keep(black_share60_lit_nc black_share60 urbanB60 unemp60 family_less_3000 pop60 school_low pro_black_county cotton_suitability cotton_share_land1964 anti_black_county rep_share_1964 urbanB60_lit unemp60_lit family_less_3000_lit pop60_lit school_low_lit cotton_suitability_lit cotton_share_land1964_lit pro_black_county_lit anti_black_county_lit rep_share_1964_lit) cells("b(star label( )  fmt(3))" " se(par label( ) fmt(3))") replace mlabels(none) collabels(none) starlevels(* 0.10 ** 0.05 *** 0.01) nostardetach label msign (--) varwidth(12) modelwidth(10) style(tex) ///
title(OLS models. Dependent Variable: Change in Black Elected Officials (1964-1980)\label{allcoefficients_controls}) ///	
prehead("\begin{table}[!h]\begin{center} \begin{threeparttable}\topcaption{@title}\renewcommand{\arraystretch}{0.5}\footnotesize \begin{tabular}{@{}l c c c c c} \toprule & (1) &(2) &(3)&(4) &(5) \\ \cmidrule(r){2-6} ") ///
prefoot(" & & & & & \\ State Trends & No & Yes & Yes & Yes & Yes \\ \midrule") ///
postfoot("\bottomrule\end{tabular} \begin{tablenotes}[flushleft] \setlength\labelsep{0pt} \scriptsize \item Robust standard errors in parenthesis clustered by judicial divisions in columns (1)-(2)-(3)-(4) and by county in column (5). ***, ** and * indicate statistical significance at the 1\%, 5\% and 10\% levels respectively. \end{tablenotes} \end{threeparttable} \end{center} \end{table}")










