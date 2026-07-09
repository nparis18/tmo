*! version 0.9.0b3-dev 2025-09-02 -- development version (program named tmo); do NOT load together with src/tmo.ado in one session (mata name collisions)

capture program drop tmo
program define tmo, eclass
    version 13

    syntax, ///
        cmd(str) x(varname) ylist(varlist) Idvar(varname) ///
        [Timevar(varname)] ///
        [LATitude(varname)] [LONgitude(varname)] [DISTTHREShold(real 0)] [miles] ///
        [THREShold(real -9)] [thresholdoff] ///
        [MISSlimit(real 0.1)] ///
        [FILEsuffix(str)] ///
        [savedyad] [load(str)] ///
        [plotq] ///
        [plothist] [plothistnbins(int 10000)] ///
        [plotse] [saveplotseest] ///
        [saveest] ///
        [scpc_cmd(str)] [scpc_uncond]
   
    ********************************
    *** RUN AND SAVE CMD OPTIONS ***
    ********************************

    // DC: Could allow alternative estimation procedures provided they 
    //   are e-class if some judicious unpacking of arguments is done? 
    qui `cmd'
    local spec `e(cmd)' 
    // DC: Adding for ereturn display
    marksample touse
    local r2a = e(r2_a)
    local r2  = e(r2_a)
    local nobs = e(N)
    local rmse = e(rmse)
    local depvar = e(depvar)

    * Assert that cmd uses supported command
    if !inlist("`spec'","regress","reghdfe","areg","ivreghdfe","ivreg2","ivregress") {
            di as error "`spec' not supported"
            exit
    }

    * Store original results and options
    local y `e(depvar)'
    scalar beta    = _b[`x']
    scalar se      = _se[`x']
    scalar N_obs   = e(N)
    scalar N_clust = e(N_clust)
    scalar df_r    = e(df_r)
    
    tempvar tmo_sample
    gen byte `tmo_sample' = e(sample)

    * Extract clusters
    local cluster `e(clustvar)'
    
    * Extract absorb vars for (iv)reghdfe adn areg
    local absorb_vars `e(absvars)' `e(absvar)'

    * Extract weights
    if "`e(wexp)'"!="" {
        local weightvar = subinstr("`e(wexp)'","=","",1)
        local weightexp [`e(wtype)'`e(wexp)']
    }
    else {
        local weightvar
        local weightexp
    }

    *** END RUN AND SAVE CMD OPTIONS ***


    *************************
    *** RUN SCPC IF GIVEN ***
    *************************
    
    if "`scpc_cmd'"!="" {
        // Install edited version of scpc that stores critical values
        // DC: Let's bundle this with tmo.  
           // This avoids issues if users have something else located in their 
        qui net install scpc, from("https://raw.githubusercontent.com/wjnkim/tmo/master/scpc_tmo") replace

        preserve
            qui keep if !missing(`longitude') & !missing(`latitude') & `tmo_sample'

            rename `latitude' s_1
            rename `longitude' s_2
            qui hashsort `idvar'
            qui `scpc_cmd'
            qui keep if e(sample)
            if "`scpc_uncond'"=="" {
                qui scpc, k(1) latlong
            }
            else {
                qui scpc, k(1) latlong uncond
            }
            scalar scpc_se=e(scpcstats)[1,2]
        
            mata: id_scpc = st_data(.,"`idvar'")
            clear

            mata: id_scpc_uniq = uniqrows(id_scpc)
            mata: Wfin_sum_vec = vec(Wfin[.,2::cols(Wfin)]*Wfin[.,2::cols(Wfin)]'):/(cols(Wfin)-1)
            mata: id_scpc_rowvec = vec(J(1,rows(Wfin),id_scpc_uniq))
            mata: id_scpc_colvec = vec(J(rows(Wfin),1,id_scpc_uniq'))
                
            mata: st_local("scpc_obsN",strofreal(rows(Wfin_sum_vec),"%50.0f"))

            gen id1=.
            gen id2=.
            gen Wfin=.
                
            qui set obs `scpc_obsN'
            
            mata: st_store(.,.,(id_scpc_rowvec,id_scpc_colvec,Wfin_sum_vec))
            
            qui keep if id1>=id2 // keep only lower triangular

            qui compress
            tempfile Wfin
            qui save `Wfin'
        restore
    }
    else {
        global scpc_cv = .
    }

    *** END RUN SCPC IF GIVEN ***



    *********************
    *** OPTION CHECKS ***
    *********************

    * Assert gtools installed
    cap which gtools
    if _rc {
        di as error "tmo requires gtools package -- please run: ssc install gtools"
        exit
    }

    * Assert that y and x are in cmd and y is the dependent variable
    local ycheck: word 2 of `cmd'
	if strpos("`cmd'","`y'")==0 | strpos("`cmd'","`x'")==0 | "`y'"!="`ycheck'" {
		di as error "cmd must contain `y' and `x' and `y' must be independent var"
		exit
	}

    * Assert that y appears only once in cmd
	if (strlen("`cmd'")-strlen(subinstr("`cmd'"," `y' ","",.)) != strlen(" `y' ")) {
		di as error "cmd contains multiple instances of `y'"
		exit
	}

    * Assert that x appears only once in cmd
	if (strlen("`cmd'")-strlen(subinstr("`cmd'","`x'","",.)) != strlen("`x'")) {
		di as error "cmd contains multiple instances of `x', please rename variables"
		exit
	}

    * Assert no duplicates in y ylist
    local done ""
    local dups ""
    foreach var in `y' `ylist' {
        confirm var `var', exact
        if strpos("`done'"," `var' ")>0 {
            local dups "`dups' `var'"
        }
        local done "`done' `var'"
    }
    if "`dups'"!="" {
        di as error "Duplicated variables in depvar/ylist: `dups'"
        exit
    }

    * Assert misslimit is between 0 and 1
    if `misslimit'<0 | `misslimit'>1 {
        di as error "misslimit() value must be between 0 and 1"
        exit
    }

    * Assert plothistnbins is positive if given
    if `plothistnbins'<=0 {
        di as error "plothistnbins() value must be positive"
        exit
    }

    * saveplotse requires plotse option
    if "`plotse'"=="" & "`saveplotse'"!="" {
        di as error "saveplotse option requires plotse option"
        exit
    }

    * Require file path for saving figures/data
    if "`filesuffix'"=="" & ("`plotse'"!="" | "`plothist'"!="" | "`savedyad'"!="" | "`saveplotseest'"!="" | "`saveest'"!="") {
        di as error "filesuffix() required for `plotse' `plothist' `savedyad' `saveplotseest' `saveest'"
    }

    * Require clustering (at least at location level) if panel
    if "`timevar'"!="" & "`cluster'"=="" {
        di as error "Clustering required for panel case (at least at `id' level)"
        exit
    }

    * Require absorb if (iv)reghdfe
    if inlist("`spec'","reghdfe","ivreghdfe", "areg") & "`absorb_vars'"=="" {
        di as error "absorb() required for `spec'"
        exit
    }
/*
    * Require xtset if xtreg
    capture xtset
    if "`spec'" == "xtreg"{
        capture xtset
        if _rc {
            di as error "ERROR: Panel no declarado. Usa `xtset´ antes de `xtreg_safe´."
            exit 198
        }
    }
*/
    * Assert idvar is unique (within timevar if panel)
	if "`idvar'"!="" & "`timevar'"=="" {
		gisid `idvar' if `tmo_sample'
	}
	if "`idvar'"!="" & "`timevar'"!="" {
		gisid `idvar' `timevar' if `tmo_sample'
	}

    * Assert longitude and latitude provided if distthreshold!=0
    if ("`longitude'"=="" | "`latitude'"=="") & `distthreshold'!=0 {
        di as error "longitude() and latitude() required for distthreshold() option"
        exit
    }

    * Assert longitude and latitude provided if scpc
    if ("`longitude'"=="" | "`latitude'"=="") & "`scpc_cmd'"!="" {
        di as error "longitude() and latitude() required for scpc_cmd() option"
        exit
    }

    * Assert scpc_cmd() if scpc_uncond option
    if "`scpc_cmd'"=="" & "`scpc_uncond'"!="" {
        di as error "scpc_cmd() required for scpc_uncond() option"
        exit
    }

    * Assert no weights if scpc
    if "`scpc_cmd'"!="" & "`weightvar'"!="" {
        di as error "weights not allowed for scpc"
        exit
    }

    * Assert distthreshold>0 if provided
    if `distthreshold'<0 {
        di as error "distthreshold() must be greater than 0"
        exit
    }

    * Assert geodist package installed if distthreshold>0
    if `distthreshold'>0 {
        cap which geodist
        if _rc {
            di as error "distthreshold() requires geodist package -- please run: ssc install geodist"
            exit
        }
    }

    * Assert custom threshold is between 0 and 1 if provided
    if `threshold'!=-9 & (`threshold'<0 | `threshold'>1) {
        di as error "Custom threshold() must be between 0 and 1"
        exit
    }

    * Store number of locations and time periods
    qui gdistinct `idvar' if `tmo_sample'
    local N = r(ndistinct)
    scalar N = `r(ndistinct)'
    mata: N = `N'
    if "`timevar'"!="" {
        qui gdistinct `timevar' if `tmo_sample'
        local T = r(ndistinct)
        scalar T = r(ndistinct)
        mata: T = `T'
    }
    else {
        local T = 1
        scalar T = 1
        mata: T = 1
    }

    *** END OPTION CHECKS ***



    *****************
    *** CLEAN CMD ***
    *****************
    
    if inlist("`spec'","reghdfe","ivreghdfe","areg") {
        * Remove any resid specified in cmd already
        local comma_start  = strpos("`cmd'", ",")
        local before_comma = substr("`cmd'", 1, `comma_start'-1)
        local after_comma  = substr("`cmd'", `comma_start'+1, .)
        local after_comma  = regexr("`after_comma'", "(res|resi|resid|residu|residua|residual|residuals)\([^)]*\)", "")
        local after_comma  = regexr("`after_comma'", "(res|resi|resid|residu|residua|residual|residuals)", "")

        * Remove saving FE options
        local after_comma = regexr("`after_comma'", "(a|ab|abs|abso|absor|absorb)\([^)]*\)", "")

        local cmd `before_comma', absorb(`absorb_vars') `after_comma'
    }

    *** END CLEAN CMD ***



    ****************************
    *** STORE ID-CLUSTER XW ****
    ****************************

    if "`cluster'"!="" & "`cluster'"!="`idvar'" {
        preserve
            qui keep if `tmo_sample'
            keep `idvar' `cluster'

            qui gduplicates drop
            cap gisid `idvar'

            if _rc {
                di as error "Only clustering by groups of locations is supported (`cluster' must be constant within `idvar')"
                exit
            }
        restore
        local clustervars `cluster'
        local clusterOff = 0
    }
    else{
        local clustervars
        local clusterOff = 1 
    }


    *** END STORE ID-CLUSTER XW ***



    *********************************
    *** STORE LOCATION DISTANCES ****
    *********************************

    if "`longitude'"!="" & "`latitude'"!="" {
        preserve
            qui keep if `tmo_sample'
            keep `idvar' `longitude' `latitude'

            qui gduplicates drop
            cap gisid `idvar'

            if _rc {
                di as error "`longitude' and/or `latitude' not constant within some `idvar'"
                exit
            }

            gen n=_n
		
            rename `latitude'  lat2
            rename `longitude' lon2
            rename `idvar' id2
            qui compress

            tempfile dist
            qui save `dist'
            
            rename id2 id1
            rename lat2 lat1
            rename lon2 lon1	
            
            qui sum n
            qui expand `r(max)'
            drop n
            
            qui hashsort id1
            by id1: gen n=_n
            
            qui merge m:1 n using `dist', assert(3) nogen
            drop n

            qui keep if id1>=id2 // keep only lower triangular
            
            qui geodist lat1 lon1 lat2 lon2, gen(dist) `miles' sphere
            qui replace dist=0 if id1==id2

            keep id1 id2 dist

            qui compress   
            tempfile dist
            qui save `dist'
        restore
    }

    *** END STORE LOCATION DISTANCES ***



    ****************************
    *** WRITE CMD FOR XTILDE ***
    ****************************

    if  inlist("`spec'","regress","reghdfe", "areg") {
        local xtildecmd = subinstr("`cmd'"," `y' "," `x' ",1)
        local xtildecmd = subinstr("`xtildecmd'"," `x' "," ",2)
        if "`spec'"=="reghdfe" {
            local xtildecmd `xtildecmd' resid
        }
    }
    if "`spec'"=="ivreghdfe" {        
        * Create xtildecmd using ivreghdfe
        local parens_start = strpos("`cmd'", "(")
        local parens_end = strpos("`cmd'", ")")
        local paren_content = substr("`cmd'", `parens_start' + 1, `parens_end' - `parens_start' - 1)
        local instr_start = strpos("`paren_content'", "=")
        local endog_part = substr("`paren_content'", 1, `instr_start'-1)
        local endog = word("`endog_part'", 1)
        if wordcount("`endog_part'") > 1 {
            di as error "Multiple endogenous regressors not supported"
            exit
        }
        if "`endog'"!="`x'" {
            di as error "Endogenous regressor `endog' is not `x'"
            exit
        }
        local instr = substr("`paren_content'", `instr_start' + 1, .)

        local remaining = subinstr("`cmd'", "ivreghdfe", "", 1)
        local remaining = subinstr("`remaining'", "`y'", "", 1)
        local comma_start = strpos("`remaining'", ",")
        local after_comma = substr("`remaining'", `comma_start'+1, .)
        local remaining = substr("`remaining'", 1, `comma_start'-1)
        local paren_start = strpos("`remaining'", "(")
        local paren_end = strpos("`remaining'", ")")
        local controls = substr("`remaining'", 1, `paren_start'-1) + " " + substr("`remaining'", `paren_end'+1, .)
        local xtildecmd1 reghdfe `endog' `instr' `controls', `after_comma'
        local xtildecmd2 reghdfe __tmo_xhat `controls', `after_comma' resid
    }
    if "`spec'"=="ivreg2" {
        * Create xtildecmd using reg
        local parens_start = strpos("`cmd'", "(")
        local parens_end = strpos("`cmd'", ")")
        local paren_content = substr("`cmd'", `parens_start' + 1, `parens_end' - `parens_start' - 1)
        local instr_start = strpos("`paren_content'", "=")
        local endog_part = substr("`paren_content'", 1, `instr_start'-1)
        local endog = word("`endog_part'", 1)
        if wordcount("`endog_part'") > 1 {
            di as error "Multiple endogenous regressors not supported"
            exit
        }
        if "`endog'"!="`x'" {
            di as error "Endogenous regressor `endog' is not `x'"
            exit
        }
        local instr = substr("`paren_content'", `instr_start' + 1, .)

        local remaining = subinstr("`cmd'", "ivreg2", "", 1)
        local remaining = subinstr("`remaining'", "`y'", "", 1)
        local comma_start = strpos("`remaining'", ",")
        local remaining = substr("`remaining'", 1, `comma_start'-1)
        local paren_start = strpos("`remaining'", "(")
        local paren_end = strpos("`remaining'", ")")
        local controls = substr("`remaining'", 1, `paren_start'-1) + " " + substr("`remaining'", `paren_end'+1, .)

        local comma_start = strpos("`cmd'", ",")
        local after_comma = substr("`cmd'", `comma_start'+1, .)
        local has_nocon = regexm("`after_comma'", "(,noc|,noco|,nocon |,nocons|,noconst|,noconsta|,noconstan|,noconstant| noc| noco| nocon| nocons| noconst| noconsta| noconstan| noconstant)")
        if `has_nocon' local nocons , nocons 
        else local nocons

        local xtildecmd1 reg `endog' `instr' `controls' `nocons'
        local xtildecmd2 reg __tmo_xhat `controls' `nocons'
    }

    if "`spec'"=="ivregress" {

        * Extract the estimation method (2sls, gmm, liml, etc.)
        local method = word("`cmd'", 2)

        * Create xtildecmd using ivreghdfe
        local parens_start = strpos("`cmd'", "(")
        local parens_end   = strpos("`cmd'", ")")
        local paren_content = substr("`cmd'", `parens_start' + 1, `parens_end' - `parens_start' - 1)

        local instr_start = strpos("`paren_content'", "=")
        local endog_part = substr("`paren_content'", 1, `instr_start'-1)
        local endog      = word("`endog_part'", 1)
        if wordcount("`endog_part'") > 1 {
            di as error "Multiple endogenous regressors not supported"
            exit
        }
        if "`endog'"!="`x'" {
            di as error "Endogenous regressor `endog' is not `x'"
            exit
        }
        local instr = substr("`paren_content'", `instr_start' + 1, .)

        local remaining = subinstr("`cmd'", "ivregress", "", 1)
        local remaining = subinstr("`remaining'", "`method'","", 1)
        local remaining = subinstr("`remaining'", "`y'", "", 1)
        local comma_start = strpos("`remaining'", ",")
        local remaining = substr("`remaining'", 1, `comma_start'-1)
        local paren_start = strpos("`remaining'", "(")
        local paren_end = strpos("`remaining'", ")")
        local controls = substr("`remaining'", 1, `paren_start'-1) + " " + substr("`remaining'", `paren_end'+1, .)

        local comma_start = strpos("`cmd'", ",")
        local after_comma = substr("`cmd'", `comma_start'+1, .)
        local has_nocon = regexm("`after_comma'", "(,noc|,noco|,nocon |,nocons|,noconst|,noconsta|,noconstan|,noconstant| noc| noco| nocon| nocons| noconst| noconsta| noconstan| noconstant)")
        if `has_nocon' local nocons , nocons 
        else local nocons

        local xtildecmd1 reg `endog' `instr' `controls' `nocons'
        local xtildecmd2 reg __tmo_xhat `controls' `nocons'
    }

    *** END WRITE CMD FOR XTILDE ***

    

    ***************
    *** RUN TMO ***
    ***************
    preserve
        if "`load'"=="" { // if dyad data already exists, can load and skip this part (programmer option)
            qui keep if `tmo_sample'

            * Estimate __tmo_xtilde
            if  inlist("`spec'","regress","reghdfe","areg") {
                qui `xtildecmd'
                predict __tmo_xtilde, resid
            }
            if inlist("`spec'","ivreghdfe","ivreg2", "ivregress") {
                qui `xtildecmd1'
                qui predict __tmo_xhat
                qui `xtildecmd2'
                qui predict __tmo_xtilde, resid
                drop __tmo_xhat
            }

            * Check __tmo_xtilde is correct
            qui reg `y' __tmo_xtilde `weightexp'
            if abs(beta-_b[__tmo_xtilde])>1e-5 {
                di as error "__tmo_xtilde is incorrect"
                exit
            }

            * Loop through auxiliary outcomes and save residuals
            cap drop __tmo_resid*

            if  inlist("`spec'","regress","ivreg2", "ivregress") {
                * Only keep cmd before comma (faster runtime)
                local comma_start = strpos("`cmd'",",")
                if `comma_start'>0 {
                    local before_comma = substr("`cmd'",1,`comma_start'-1)
                    local after_comma = substr("`cmd'", `comma_start'+1, .)

                    * Check whether there is nocons option and include if so
                    local has_nocon = regexm("`after_comma'", "(,noc|,noco|,nocon |,nocons|,noconst|,noconsta|,noconstan|,noconstant| noc| noco| nocon| nocons| noconst| noconsta| noconstan| noconstant)")
                    if `has_nocon' {
                        local nocons , nocons 
                    }
                    else local nocons
                }
                else {
                    local before_comma `cmd'
                    local nocons
                }
                
                local cmd_toloop `before_comma' `nocons'

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    qui `cmd_inloop'
                    qui predict __tmo_resid`ynum', resid
                    local ++ynum
                }
            }
            
            if "`spec'"=="reghdfe" {
                * Remove any clustering (faster runtime)
                //Note DC: build in catch for spaces between option name and parentheses
                local cmd_toloop `cmd'
                while regexm("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)") {
                    local cmd_toloop = regexr("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)", "")
                }

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    qui `cmd_inloop' resid
                    qui predict __tmo_resid`ynum', resid
                    local ++ynum
                }
            }
            if "`spec'"=="areg" {
                * Remove any clustering (faster runtime)
                local cmd_toloop `cmd'
                while regexm("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)") {
                    local cmd_toloop = regexr("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)", "")
                }

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    qui `cmd_inloop'
                    qui predict __tmo_resid`ynum', resid
                    local ++ynum
                }
            }
            if "`spec'"=="ivreghdfe" {
                * Remove any clustering (faster runtime)
                local cmd_toloop `cmd'
                while regexm("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)") {
                    local cmd_toloop = regexr("`cmd_toloop'", "(cl|clu|clus|clust|cluste|cluster|vce)\([^)]*\)", "")
                }

                * Specify FEs to save
                local absorb_start = regexm("`cmd'", "(a|ab|abs|abso|absor|absorb)\(([^)]+)\)")
                local absorb_vars "`=regexs(2)'"
                local absorb_vars_savefe
                local absorb_vars_fesum
                local fe=1
                foreach fe_var in `absorb_vars' {
                    local absorb_vars_savefe `absorb_vars_savefe' __tmo_fe`fe'=`fe_var'
                    local absorb_vars_fesum `absorb_vars_fesum' + __tmo_fe`fe'
                    local ++fe
                }
                local cmd_toloop = regexr("`cmd_toloop'", "(a|ab|abs|abso|absor|absorb)\([^)]*\)", "")
                local cmd_toloop `cmd_toloop' absorb(`absorb_vars_savefe')

                local ynum=1
                foreach aux_y in `y' `ylist' {
                    local cmd_inloop = subinstr("`cmd_toloop'","`y'","`aux_y'",1)
                    cap drop __tmo_fe*
                    qui `cmd_inloop'
                    cap drop __tmo_xb 
                    qui predict __tmo_xb
                    qui gen __tmo_resid`ynum' = `aux_y' - (__tmo_xb`absorb_vars_fesum')
                    local ++ynum
                }
            }
            scalar D = `ynum'-1

            * Calculate correlation in residuals and contribution to variance
            keep `idvar' `timevar' `weightvar' __tmo_xtilde __tmo_resid* `clustervars'
            qui hashsort `idvar' `timevar'

            * Store data in Mata
            mata: id = st_data(.,"`idvar'")
            mata: xtilde = st_data(.,"__tmo_xtilde")
            
            if "`weightvar'"!="" {
                mata: wgt = st_data(.,"`weightvar'")
            }
            else {
                mata: wgt = J(rows(xtilde),1,1)
            }
            mata: xtilde_wgt = xtilde:*wgt
            mata: Res1 = st_data(.,"__tmo_resid1")

            mata: CovEpsVec = J(0,0,.)
            mata: id_widerowvec = J(0,0,.)
            mata: id_widecolvec = J(0,0,.)
            mata: res1_widerowvec = J(0,0,.)
            mata: res1_widecolvec = J(0,0,.)
            mata: xtilde_widerowvec = J(0,0,.)
            mata: xtilde_widecolvec = J(0,0,.)

            if "`timevar'"=="" { // Cross-sectional case
                * Compute correlation in outcomes
                //DC: Minor bottleneck here, but not too bad
                mata: Res = st_data(.,"__tmo_resid*")

                if(`clusterOff') {
                    clear
                    mata: corr_resid(id, Res, Res1, `misslimit', CovEpsVec, xtilde_wgt, id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, N)
                }
                else {
                    mata: id_cluster = st_data(.,"`clustervars'")
                    clear
                    mata: same_cl = J(0,0,.)
                    mata: corr_resid(id, Res, Res1, `misslimit', CovEpsVec, xtilde_wgt, id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, N, id_cluster, same_cl)
                }
                mata: st_local("obsN",strofreal(rows(CovEpsVec),"%50.0f"))
                gen id1=.
                gen id2=.
                gen corr=.

                qui set obs `obsN'

                mata: st_store(.,.,(id_widerowvec,id_widecolvec,CovEpsVec))

                qui gen corr_fisher = 0.5 * ln((1+corr)/(1-corr))
                cap drop offdiag
                qui gen byte offdiag = !missing(corr) & id1!=id2 & abs(corr)!=1

                if `threshold' == -9{
                    qui gstats sum corr_fisher if offdiag
                    scalar sd = (r(p75)-r(p25))/(invnormal(0.75)-invnormal(0.25))
                    scalar df = 1/(sd^2)
                    cap drop corr_fisher_abs 
                    qui gen corr_fisher_abs = abs(corr_fisher)
                    qui replace corr_fisher_abs = . if !offdiag // will be at end
                    qui hashsort corr_fisher_abs

                    cap drop cdf_emp_abs_1min cdf_iqr_abs_1min q_iqr_abs
                    qui count if offdiag
                    qui gen cdf_emp_abs_1min = 1 - (_n/`r(N)') if offdiag
                    qui gen cdf_iqr_abs_1min = 1 - 2*(normal(corr_fisher_abs/sd)-0.5) if offdiag
                    qui gen q_iqr_abs = cdf_emp_abs_1min - 2*cdf_iqr_abs_1min if offdiag

                    qui sum q_iqr_abs if offdiag
                    qui sum corr_fisher_abs if abs(q_iqr_abs-`r(max)')<=1e-10 & offdiag
                    scalar fthres = `r(min)'
                    scalar thres = tanh(fthres)
                    local thres = thres
                    keep id1 id2 corr corr_fisher corr_fisher_abs q_iqr_abs offdiag
                }
                else{
                    local thres = `threshold'
                    keep id1 id2 corr corr_fisher offdiag
                }
               * Compute contribution to SE for each pair of locations
                mata: ResXtildeVec = J(0,0,.)

                //DC: THIS IS THE BOTTLENECK
                if(`clusterOff') {
                    mata: sandwich_crosssec(id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, CovEpsVec, ResXtildeVec, `thres')
                    local posSW = 3
                }
                else{
                    mata: sandwich_crosssec(id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, CovEpsVec, ResXtildeVec, `thres', same_cl)
                    local posSW = 4
                }
            }
            else { // Panel case
                * Store data in Mata
                mata: t = st_data(.,"`timevar'")

                ** Reshape Res
                * To make sure each location is shown for all time periods
                qui  tsset `idvar' `timevar'
                tsfill, full

                * Drop time periods that are all missing for resids of main outcome 
                qui gegen __tmo_resid1_missing = sum(!missing(__tmo_resid1)), by(`timevar')
                qui drop if __tmo_resid1_missing==0
                drop __tmo_resid1_missing
                
                * Input NT x D matrix of residuals to Mata
                mata: Res = st_data(.,"__tmo_resid*")
                mata: Dn = cols(Res)
                mata: idfull = st_data(.,"`idvar'")
 
                * Reshape matrix of residuals to N x TD
                mata: Res = rowshape(Res,N)
                mata: Res1= rowshape(Res1,N)
                mata: id_wide = rowshape(idfull,N)
                mata: id_wide = id_wide[,1]
                mata: assert(rows(uniqrows(id_wide))==N)
                mata: assert(cols(Res)==Dn*T)

                if(`clusterOff') {
                    clear
                    mata: corr_resid(id_wide, Res, Res1, `misslimit', CovEpsVec, xtilde_wgt, id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, N)
                    local posSW = 3
                }
                else {
                    mata: id_cluster = st_data(.,"`clustervars'")
                    clear
                    mata: same_cl = J(0,0,.)
                    mata: corr_resid(id_wide, Res, Res1, `misslimit', CovEpsVec, xtilde_wgt, id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, N, id_cluster, same_cl)
                    local posSW = 4
                }
                * Compute correlation in outcomes
                mata: st_local("obsN",strofreal(rows(CovEpsVec),"%50.0f"))
                gen id1=.
                gen id2=.
                gen corr=.

                qui set obs `obsN'

                mata: st_store(.,.,(id_widerowvec,id_widecolvec,CovEpsVec))

                qui gen corr_fisher = 0.5 * ln((1+corr)/(1-corr))
                cap drop offdiag
                qui gen byte offdiag = !missing(corr) & id1!=id2 & abs(corr)!=1

                qui gstats sum corr_fisher if offdiag
                scalar sd = (r(p75)-r(p25))/(invnormal(0.75)-invnormal(0.25))
                scalar df = 1/(sd^2)
                cap drop corr_fisher_abs 
                qui gen corr_fisher_abs = abs(corr_fisher)
                qui replace corr_fisher_abs = . if !offdiag // will be at end
                qui hashsort corr_fisher_abs

                cap drop cdf_emp_abs_1min cdf_iqr_abs_1min q_iqr_abs
                qui count if offdiag
                qui gen cdf_emp_abs_1min = 1 - (_n/`r(N)') if offdiag
                qui gen cdf_iqr_abs_1min = 1 - 2*(normal(corr_fisher_abs/sd)-0.5) if offdiag
                qui gen q_iqr_abs = cdf_emp_abs_1min - 2*cdf_iqr_abs_1min if offdiag

                qui sum q_iqr_abs if offdiag
                qui sum corr_fisher_abs if abs(q_iqr_abs-`r(max)')<=1e-10 & offdiag
                scalar fthres = `r(min)'
                           
                if `threshold' == -9 {
                    scalar thres = tanh(fthres)
                    local thres = thres
                }
                else local thres = `threshold'
                keep id1 id2 corr corr_fisher corr_fisher_abs q_iqr_abs offdiag

                * Compute contribution to SE for each pair of locations
                sandwich_panel, clusterOff(`clusterOff') nloc(`N') ntime(`T') th(`thres') noi
            }
            * Normalize contribution to variance
            mata: resxtildenorm(ResXtildeVec, xtilde_wgt, `posSW')

            qui compress
            tempfile corr
            qui save `corr'
            local corr corrpath(`corr')

            clear

            * Calculate part of contribution to scpc SE if specified
            if "`scpc_cmd'"!="" {
                mata: y_scpc = (sqrt(N)/(xtilde'*xtilde)):*Res1:*xtilde	
                mata: st_local("scpc_obsN",strofreal(rows(id),"%50.0f"))
                
                clear
                gen id1=.
                gen y_scpc1=.
                
                qui set obs `scpc_obsN'
                
                mata: st_store(.,.,(id,y_scpc))
                
                cap gisid id1
                if _rc gcollapse (sum) y_scpc1, by(id1)

                qui compress
                tempfile scpc1
                qui save `scpc1'
                
                rename id1 id2
                rename y_scpc1 y_scpc2
                
                tempfile scpc2
                qui save `scpc2'
            }
            if "`filesuffix'"!="" {
                local savepath filesuffix(`filesuffix') 
            }
            else {
                local savepath
            }
            if "`longitude'"!="" & "`latitude'"!="" {
                local distp distpath(`dist')
            }
            else {
                local distp
            }
            if "`scpc_cmd'"!="" {
                local scpcw scpcwfin(`Wfin')
                local scpcy1 scpcy1path(`scpc1')
                local scpcy2 scpcy2path(`scpc2')
            }
            else {
                local scpcw
                local scpcy1
                local scpcy2
            }
            load_data, nloc(`N') clusterOff(`clusterOff') `distp' `savedyad' `savepath' `scpcw' `scpcy1' `scpcy2' `corr'
        }
        else {
            if "`filesuffix'"!="" {
                local savepath filesuffix(`filesuffix') 
            }
            else {
                local savepath
            }

            use "`load'_est.dta", clear
            scalar D = N_outcomes[1]

            use "`load'_dyad.dta", clear
        }
        * Compute TMO SE
        est_tmo_se, dist_cutoff(`distthreshold') `thresholdoff'

        if ("`plotq'"!="" | "`plothist'"!="" ) {
            cap drop pdf_iqr
            qui gen pdf_iqr = normalden(corr_fisher,sd) if offdiag
            plots_dfqt, `plotq' `plothist' nbins(`plothistnbins') nloc(`N') `savepath'
        }

        * Plot TMO SE over threshold 
        if "`plotse'"!="" {
            tmo_over_thres, `savepath' `saveplotseest' noi dist_cutoff(`distthreshold')
        }   
    restore

        **********************
        *** OUTPUT RESULTS ***
        **********************
              
            mat tmo_results = J(1,6,.)
            mat rownames tmo_results = "`x'"
            mat colnames tmo_results = "Coef" "TMO SE" "t" "P>|t|" "95% Conf" "Interval"
            mat tmo_results[1,1] = beta
            mat tmo_results[1,2] = tmo_se
            mat tmo_results[1,3] = beta/tmo_se
            mat tmo_results[1,4] = 2*ttail(df_r, abs(beta/tmo_se))
            scalar lb = beta - invttail(df_r,0.025)*tmo_se
            scalar ub = beta + invttail(df_r,0.075)*tmo_se
            mat tmo_results[1,5] = lb
            mat tmo_results[1,6] = ub
            //matlist tmo_results, border(all) cspec(o2& %20s | %9.3f o2 & %9.3f o2 & %6.2f o2 & %4.3f o2 & %9.3f o2 & %9.3f o2 &) rspec(&-&)

            mat tmo_details = J(5,1,.)
            mat rowname tmo_details = "Optimal threshold" "% of off-diag in SE est." "% >= threshold (excl. clusters/Conley)" "# outcomes" "Degrees of freedom" 
            mat tmo_details[1,1] = thres
            mat tmo_details[2,1] = offdP*100
            mat tmo_details[3,1] = offdPnocl*100
            mat tmo_details[4,1] = D
            mat tmo_details[5,1] = df
            
            matrix b = beta
            matrix V = tmo_se^2

            matrix colnames b = `x'
            matrix colnames V = `x'
            matrix rownames V = `x'
            
            local df = df_r
            ereturn post b V, dof(`df') obs(`nobs') esample(`touse') depname(`depvar')

            ereturn scalar r2 = `r2'
            ereturn scalar r2_a = `r2a'
            ereturn scalar rmse = `rmse'
            _coef_table_header, nomodeltest title(Linear Regression with TMO)
            dis ""
            _coef_table
            matlist tmo_details, cspec(& %-68s & %9.3f &) rspec(- & & & & & -) coleqonly noblank

            //ereturn display, level(`level')
            
            ereturn scalar beta = beta
            ereturn scalar orig_se = se
            ereturn scalar tmo_se = tmo_se
            ereturn scalar lb = lb
            ereturn scalar ub = ub
            ereturn scalar threshold = thres
            ereturn scalar pct_ge_thres = offdP*100
            ereturn scalar pct_ge_thres_nocl = offdPnocl*100
            ereturn scalar T = T
            ereturn scalar N_loc = N
            ereturn scalar N_clust = N_clust
            ereturn scalar N_outcomes = D
            ereturn scalar N = N_obs
            ereturn scalar dof = df
            ereturn scalar finite_sample_dof = dof_adj
            ereturn scalar df_r = df_r
            ereturn scalar scpc_cv = ${scpc_cv}

        if "`saveest'"!="" {
            tmo_save, `savepath'
        }
    end


***********************
*** HELPER COMMANDS ***
***********************

* Function for computing correlation of residuals
cap mata: mata drop corr_resid()
mata
	void corr_resid(real matrix id,
                    real matrix Res, 
                    real matrix Res1, 
                    real scalar misslimit, 
                    real matrix CovEpsVec,
                    real matrix xtilde, 
                    real matrix id_widerowvec, 
                    real matrix id_widecolvec, 
                    real matrix res1_widerowvec, 
                    real matrix res1_widecolvec, 
                    real matrix xtilde_widerowvec,
                    real matrix xtilde_widecolvec,
                    real scalar N,
                    |real matrix cluster, 
                    real matrix same_cl)
	{
		real matrix Res_ms, Res_ms_lethres, Res_ms_ind, Res_ms_ind_tr, Res_no_ms, Denom, ResSum_DinBoth, ResMean_DinBoth, ResMeanProd_DinBoth, CovEps, Sum_ResSq_DinBoth, ResSD_DinBoth, DenomCorr, DenomVec, LowTriInd, rowpos, colpos, idx
		real scalar demean

        // Drop outcomes that are missing for more than `misslimit'
        Res_ms = colsum(Res:==.)
        Res_ms_lethres = (Res_ms:/rows(Res)):<=misslimit
        Res = select(Res,Res_ms_lethres)

        // standardize residuals
        Res = Res:-J(rows(Res),1,colsum(Res):/colsum(Res:!=.))
        Res = Res:/J(rows(Res),1,(colsum(Res:^2):/colsum(Res:!=.)):^0.5) // Studentize residuals

        assert (sum(abs(colsum(Res)):<1e-5)==cols(Res))
        assert (sum(abs((colsum(Res:^2):/colsum(Res:!=.)):^0.5 :- 1):<1e-5)==cols(Res))

        // covariance of residuals -> correlations
        Res_ms_ind = Res:!=. // keep track of missing (0==missing)
        Res_ms_ind_tr = Res_ms_ind'
        Res_no_ms = editmissing(Res,0) // Res with missing replaced with 0
        Denom = Res_ms_ind*Res_ms_ind_tr // number of both nonmissing        

        ResSum_DinBoth = Res_no_ms * Res_ms_ind_tr
        ResMean_DinBoth = ResSum_DinBoth :/ Denom
        ResMeanProd_DinBoth = ResMean_DinBoth :* ResMean_DinBoth'
        
        demean=1 // set to 1 to make residuals mean 0 within each location
        if (demean==1) {
            CovEps = ((Res_no_ms*Res_no_ms'):/Denom) - (ResMean_DinBoth:*(ResSum_DinBoth':/Denom)) - (ResMean_DinBoth':*(ResSum_DinBoth:/Denom)) + ResMeanProd_DinBoth
        }
        else {
            CovEps = ((Res_no_ms*Res_no_ms'):/Denom)
        }
        
        Sum_ResSq_DinBoth = (Res_no_ms:^2) * Res_ms_ind_tr
        ResSD_DinBoth = ((Sum_ResSq_DinBoth:/Denom) + (-2:*ResMean_DinBoth:*(ResSum_DinBoth:/Denom)) + ResMean_DinBoth:^2) :^ 0.5
        
        DenomCorr = ResSD_DinBoth :* ResSD_DinBoth'
        
        if (demean==1) {
            CovEps = CovEps:/DenomCorr
        }
        else {
            CovEps = CovEps :/ ((diagonal(CovEps):^0.5) * (diagonal(CovEps):^0.5)')
        }

        // vectorize off-diag covariance of residuals
        CovEpsVec = vec(CovEps)
        LowTriInd = J(rows(CovEps),cols(CovEps),1)
        _lowertriangle(LowTriInd,1)
        LowTriInd = vec(LowTriInd)
        CovEpsVec = select(CovEpsVec,LowTriInd)

        assert (rows(CovEpsVec)==rows(Res)*(rows(Res)+1)/2)
        
        // vectorize Denom
        DenomVec = vec(Denom)
        DenomVec = select(DenomVec,LowTriInd)

        // store id vectors
        idx     = selectindex(LowTriInd)
        rowpos  = mod(idx :- 1, N) :+ 1
        colpos  = floor((idx :- 1) :/ N) :+ 1

        id_widerowvec    = id[rowpos, .]
        id_widecolvec    = id[colpos, .]
        res1_widerowvec  = Res1[rowpos, .]
        res1_widecolvec  = Res1[colpos, .]
        xtilde_widerowvec  = (rowshape(xtilde,N))[rowpos, .]
        xtilde_widecolvec  = (rowshape(xtilde,N))[colpos, .]
        if(args()==15){
            same_cl  = ((rowshape(cluster,N))[rowpos, .])[.,1]:== ((rowshape(cluster,N))[colpos, .])[.,1]
        }
	}
end

* Function for computing contribution to SE for each pair of locations in cross-sectional case
cap mata: mata drop sandwich_crosssec()
mata
    void sandwich_crosssec(real matrix id_widerowvec, 
                           real matrix id_widecolvec,
                           real matrix res1_widerowvec, 
                           real matrix res1_widecolvec, 
                           real matrix xtilde_widerowvec, 
                           real matrix xtilde_widecolvec, 
                           real matrix CovEpsVec,
                           real matrix ResXtildeVec,
                           real scalar thres,
                         | real matrix same_cl)
	{   
        real matrix IdRes1Xtilde, res_xtide1, keep_th

        if (args() == 9) {
            keep_th = abs(CovEpsVec) :>= thres
            IdRes1Xtilde  = select((id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec), keep_th)
            ResXtildeVec  = (IdRes1Xtilde[.,3]:* IdRes1Xtilde[.,4]):*(IdRes1Xtilde[.,5]:* IdRes1Xtilde[.,6])
            ResXtildeVec  = (IdRes1Xtilde[.,1..2], ResXtildeVec)
            ResXtildeVec
        }
        else{
            real keep_cl
            keep_th = (abs(CovEpsVec) :>= thres) :| same_cl
            IdRes1Xtilde  = select((id_widerowvec, id_widecolvec, res1_widerowvec, res1_widecolvec, xtilde_widerowvec, xtilde_widecolvec, same_cl), keep_th)
            ResXtildeVec  = (IdRes1Xtilde[.,3]:* IdRes1Xtilde[.,4]):*(IdRes1Xtilde[.,5]:* IdRes1Xtilde[.,6])
            ResXtildeVec  = (IdRes1Xtilde[.,1], IdRes1Xtilde[.,2], IdRes1Xtilde[.,7], ResXtildeVec)
        }
    }
end

* Functions for computing contribution to SE for each pair of locations in panel case
cap program drop sandwich_panel
program define sandwich_panel
    syntax , clusterOff(int) nloc(int) ntime(int) th(real) [NOIsily]

    preserve
    if(`clusterOff') {
        mata: keep_th = abs(CovEpsVec) :>= `th'
        local cl=2
        local colsR1=2*`ntime' + `cl'
        local colsXtF=`colsR1' + 1
        local colsXtL=`colsR1' + 2*`ntime' 
        mata: Mout = J(0,3,.)
        mata: allSample = select((res1_widerowvec, res1_widecolvec, id_widerowvec, id_widecolvec, xtilde_widerowvec, xtilde_widecolvec), keep_th)
    }
    else {
        mata: keep_th = abs(CovEpsVec) :>= `th' :| same_cl
        local cl=3
        local colsR1=2*`ntime' + `cl'
        local colsXtF=`colsR1' + 1
        local colsXtL=`colsR1' + 2*`ntime' 
        mata: Mout = J(0,4,.)
        mata: allSample = select((res1_widerowvec, res1_widecolvec, id_widerowvec, id_widecolvec, same_cl, xtilde_widerowvec, xtilde_widecolvec), keep_th)
    }

    mata: res1all   = allSample[.,1..`colsR1']
    mata: xtall     = allSample[.,`colsXtF'..`colsXtL']
    mata: P         = rows(res1all)
    mata: st_local("P", strofreal(P, "%12.0f"))

    // We replicate the cross-section case in terms of RAM usage
    local chunk = `P'/`ntime'

    forvalues start = 1(`chunk')`P' {
        local end = min(`start' + `chunk' - 1, `P')

        if "`noisily'"!="" {
            local donepct = string(`start'*100/`P', "%5.2f")
            di "Computed `donepct'% of sandwich"
        }
        mata: sandwich_panel_loop(`start',`end',`ntime', ///
                                   res1all, xtall, Mout, `cl')
    }
    mata: ResXtildeVec = Mout
    mata: st_numscalar("num_pairs", rows(ResXtildeVec))
    restore
end

cap mata: mata drop sandwich_panel_loop()
mata
void sandwich_panel_loop(real scalar start_i,
                         real scalar end_i,
                         real scalar T,
                         real matrix Mres1,
                         real matrix Mxt,
                         real matrix Mout,
                         real scalar cl)
{
    real colvector Mout1
    real matrix sub, usub, xsub, x_u, block, rowSumx_u
    real scalar i, pairT, pair1, pair2

    pairT= 2*T
    pair1=pairT + 1
    pair2=pairT + cl

    sub  = (Mres1)[start_i..end_i, .]
    usub = sub[., 1..pairT]
    xsub = Mxt[start_i..end_i , .]
    x_u  = usub:*xsub
    rowSumx_u = rowsum(x_u[., T+1..pairT])
    block = J(rows(x_u), T, .)

    for (i = 1; i <= T; i++) {
        block[., i] = (x_u[., i] :* rowSumx_u)
    }
    Mout1 = rowsum(block)

    Mout = Mout \ ( sub[., pair1..pair2], Mout1)
}
end

* Function to normalize contribution to SE by 2*[X'X]^-1 and multiply off-diag by 2 (since only lower triangular)
cap mata: mata drop resxtildenorm()
mata
    void resxtildenorm(real matrix ResXtildeVec, 
                       real matrix xtilde_wgt,
                       real scalar posSW)
	{
		real matrix denom, offdiag
        real colvector v

        v = ResXtildeVec[.,posSW]
        denom = 1/(quadcross(xtilde_wgt,xtilde_wgt))^2
		offdiag = ResXtildeVec[.,1]:!=ResXtildeVec[.,2]
        ResXtildeVec[.,posSW] = (v :* denom) :* (1 :+ offdiag)
    }
end

* Function to load Mata data into Stata
cap program drop load_data
program define load_data,
    syntax, nloc(int) clusterOff(int) [distpath(str)] [savedyad] [filesuffix(str)] [scpcwfin(str)] [scpcy1path(str)] [scpcy2path(str)] [corrpath(str)]

    clear
    mata: st_local("obsN",strofreal(rows(ResXtildeVec),"%50.0f"))

    gen id1=.
    gen id2=.
    qui set obs `obsN'

    qui if `clusterOff'==0 {
        gen same_cl=.
        replace same_cl=0 if same_cl==.
    }
    qui gen xxresxx=.
    mata: st_store(.,.,(ResXtildeVec))
    qui merge 1:1 id1 id2 using `corrpath', assert(2 3) nogen

    qui compress

    if _N!=(`nloc'^2+`nloc')/2 {
        di as error "Number of location pairs is incorrect"
        exit
    }

    if "`distpath'"!="" {
        qui merge 1:1 id1 id2 using `distpath', assert(3) nogen
    }

    if "`scpcwfin'"!="" {
        qui merge 1:1 id1 id2 using `scpcwfin', assert(1 3) nogen
        qui merge m:1 id1 using `scpcy1path', assert(1 3) nogen
        qui merge m:1 id2 using `scpcy2path', assert(1 3) nogen
        qui gen ryyr = (Wfin*y_scpc1*y_scpc2)*(1 + (id1!=id2))
    }

    * Fisher-transform correlation
    qui count if abs(corr)>1 & !missing(corr)
    if `r(N)'>0 {
        di as error "Correlations <-1 or >1 exist"
        exit
    }
    if "`savedyad'"!="" {
        save "`filesuffix'_dyad.dta", replace
    }
end


* Function to estimate TMO SE
cap program drop est_tmo_se
program define est_tmo_se,
    syntax, [dist_cutoff(real 0)] [custom_thres(real -9)] [thresholdoff] [`scpc_uncond']

    if `custom_thres'!=-9 {
        scalar thres = `custom_thres'
        scalar fthres = 0.5 * ln((1+thres)/(1-thres))
        scalar df = .
        scalar sd = .
    }

    if "`thresholdoff'"!="" {
        scalar thres = .
        scalar fthres = .
        scalar df = .
        scalar sd = .
    }
    cap confirm var same_cl
    if !_rc {
        local orig_cond (same_cl)
    }
    else {
        local orig_cond (id1==id2)
    }
    cap confirm var same_cl
    if !_rc {
        if `dist_cutoff'>0 {
            local keep_cond (same_cl | dist<=`dist_cutoff')
        }
        else {
            local keep_cond (same_cl)
        }
    }
    else {
        if `dist_cutoff'>0 {
            local keep_cond (id1==id2 | dist<=`dist_cutoff')
        }
        else {
            local keep_cond (id1==id2)
        }
    }

    * For SCPC option
    cap confirm var Wfin
    if !_rc {
        local scpc = 1
        * Check SCPC SE (might not equal if uncond option or if missing some locations due to missing coordinates)
        qui count if missing(ryyr)
        if r(N)==0 & "`scpc_uncond'"=="" {
            qui sum ryyr
            if abs(scpc_se - sqrt(r(sum)))>1e-5 {
                di as error "SCPC SE does not match"
                exit
            }
        }
    }
    else {
        local scpc = 0
    }

    * Back out finite sample degrees of freedom adjustment from original SE
    qui sum xxresxx if `orig_cond'
    scalar dof_adj = se^2/r(sum)

    * TMO SE
    qui sum xxresxx
    local xxresxx_sum = r(sum)
    if `scpc'==1 {
        qui sum ryyr if (abs(corr)<thres | missing(corr)) & !`keep_cond'
        local ryyr_sum = r(sum)
    }
    else {
        local ryyr_sum = 0
    }
	scalar tmo_se = sqrt((`xxresxx_sum'+`ryyr_sum')*dof_adj)
    * Store no. off-diag
    qui count if id1!=id2
    scalar offdN = r(N)
    qui count if !`orig_cond'
    scalar offdNnocl = r(N)

    * Proportion off-diag included in SE calculation (including both clustering, Conley, and TMO)
	qui count if (((abs(corr)>=thres) & !missing(corr)) | `keep_cond') & id1!=id2
	scalar offdP = r(N)/offdN

    * Proportion off-diag over threshold outside clusters or Conley
	qui count if (abs(corr)>=thres) & !missing(corr) & !`keep_cond'
	scalar offdPnocl = r(N)/offdNnocl
end

* Function to plot TMO SE over thresholds
cap program drop tmo_over_thres
program define tmo_over_thres,
    syntax, filesuffix(str) [saveplotseest] [NOIsily] [dist_cutoff(real 0)]

    cap confirm var same_cl
    if !_rc {
        if `dist_cutoff'>0 {
            local keep_cond (same_cl | dist<=`dist_cutoff')
            local subtitle2 `" "(excluding within-cluster and -`dist_cutoff' correlations)" "'
        }
        else {
            local keep_cond (same_cl)
            local subtitle2 `" "(excluding within-cluster correlations)" "'
        }
    }
    else {
        if `dist_cutoff'>0 {
            local keep_cond (id1==id2 | dist<=`dist_cutoff')
            local subtitle2 `" "(excluding within-`dist_cutoff' correlations)" "'
        }
        else {
            local keep_cond (id1==id2)
            local subtitle2
        }
    }

    cap confirm var Wfin
    if !_rc {
        local scpc = 1
    }
    else {
        local scpc = 0 
    }

    * Initialize matrix to store results
    mat tmo_over_thres = J(101,4,.)
	mat colnames tmo_over_thres = delta tmo_se offdP offdPnocl 

    local row=1
    forv thr=0(0.01)1.01 {
        if "`noisily'"!="" {
            if mod(`row',20)==1 {
                local thr_str = string(`thr',"%03.2f")
                di "Calculating TMO SE over threshold at `thr_str'"
            }
        }

        mat tmo_over_thres[`row',1]=`thr'
		qui sum xxresxx if ((abs(corr)>=`thr') & !missing(corr)) | `keep_cond'
        local xxresxx_sum = r(sum)
        if `scpc'==1 {
            qui sum ryyr if (abs(corr)<thres | missing(corr)) & !`keep_cond'
            local ryyr_sum = r(sum)
        }
        else {
            local ryyr_sum = 0
        }    
		mat tmo_over_thres[`row',2]=sqrt((`xxresxx_sum'+`ryyr_sum')*dof_adj)
		qui count if (abs(corr)>=`thr') & !missing(corr) & !(id1==id2)
		mat tmo_over_thres[`row',3]=r(N)/offdN
		qui count if (abs(corr)>=`thr') & !missing(corr) & !(`keep_cond')
		mat tmo_over_thres[`row',4]=r(N)/offdNnocl

        local ++row
    }

    * Plot TMO SE over threshold
    clear
    qui svmat2 tmo_over_thres, names(col)
    qui compress

    if "`saveplotseest'"!="" {
        save "`filesuffix'_tmo_se_over_thres.dta", replace
    }

    qui gen tmo_orig_se_ratio = tmo_se/se
    qui replace tmo_orig_se_ratio = 0 if missing(tmo_orig_se_ratio)

    local thres = thres
    twoway 	(line tmo_orig_se_ratio delta, lcolor(black) lwidth(medthick)) ///
			, ///
			graphregion(color(white)) ///
			xlabel(0(0.1)1.0, format("%02.1f") grid gmin gmax) ///
			ylab(, angle(horizontal) gmin gmax) ///
			xtitle("Threshold {it:{&delta}}") ///
			subtitle("Ratio of TMO standard error to original", pos(11)) ytitle("") ///
			legend(off) ///
			yline(1, lcolor(gray%50) lwidth(thick)) ///
			xline(`thres', lcolor(red) lwidth(medthick)) ///
			xsize(16) ysize(9)
	graph export "`filesuffix'_se_ratio_over_thres.pdf", as(pdf) replace

    twoway 	(line offdPnocl delta if offdPnocl<=0.1, lcolor(black) lwidth(medthick)) ///
			, ///
			graphregion(color(white)) ///
			ylab(0(0.01)0.1, format("%03.2f") angle(horizontal) gmin gmax) ///
			xlabel(0(0.1)1.0, format("%02.1f") grid gmin gmax) ///
			xtitle("Threshold {it:{&delta}}") ///
			subtitle("Proportion of correlations {&ge} threshold" `subtitle2', pos(11)) ytitle("") ///
			legend(off) ///
			xline(`thres', lcolor(red) lwidth(medthick)) ///
			xsize(16) ysize(9)
	graph export "`filesuffix'_prop_above_thres.pdf", as(pdf) replace
end

cap program drop plots_dfqt
program define plots_dfqt
    syntax, [plotq] [plothist] nbins(int) nloc(int) [filesuffix(str)]

    if "`plotq'"!="" {        
        qui sum corr_fisher_abs if q_iqr_abs>=-0.002 & offdiag
        local xstart=floor(`r(min)'*10)/10
        local fthres=fthres

        twoway ///
                (line q_iqr_abs corr_fisher_abs if q_iqr_abs>=-0.002 & corr_fisher_abs<=1, lcolor(blue) lwidth(medthick) sort(corr_fisher_abs)), ///
                graphregion(color(white)) ///
                yline(0, lcolor(gray)) ///
                xtitle("Threshold {it:{&delta}}") ytitle("") ///
                ylab(, angle(horizontal) format("%04.3f")) ///
                xlab(`xstart'(0.1)1, grid gmin gmax format("%02.1f")) ///
                xline(`fthres', lcolor(red)) ///
                xsize(16) ysize(9)
        graph export "`filesuffix'_qt.pdf", as(pdf) replace
    }

    if "`plothist'"!="" {
        qui sum corr_fisher if offdiag
        local binwidth = (`r(max)'-(`r(min)'))/`nbins'
        qui gen bin_corr = floor((corr_fisher-`r(min)')/`binwidth') if offdiag
        qui gen bin_cent = bin_corr*`binwidth' + `binwidth'/2 + `r(min)'
        qui sum bin_corr if offdiag
        assert abs(`r(max)'-`nbins')<=1
        qui hashsort bin_corr 
        
        qui by bin_corr: gen binN=_N if offdiag
        qui gen bin_dens = binN/(_N*`binwidth') if offdiag
        qui gegen byte bin_tag = tag(bin_corr) if offdiag
        
        local thres_str = string(thres,"%03.2f")
        local fthres_str = string(fthres,"%03.2f")
        local df_str = string(df,"%05.2f")
        local fthres = fthres
        
        twoway 	(bar bin_dens bin_cent if bin_tag==1 & corr_fisher>=-1 & corr_fisher<=1, base(0) barwidth(`binwidth') color(midgreen%30)) ///
                    (line pdf_iqr corr_fisher if corr_fisher>=-1 & corr_fisher<=1, sort(corr_fisher) lcolor(blue%90)) ///
                    , graphregion(color(white)) ///
                    xtitle("Fisher transformed correlation") ///
                    ytitle("Density") ///
                    ylab(, format(%02.1f) angle(horizontal)) ///
                    xlab(-1(0.2)1, format(%02.1f)) ///
                    xline(0, lcolor(gray) lpattern(longdash)) ///
                    xline(`fthres', lcolor(red)) ///
                    xline(-`fthres', lcolor(red)) ///
                    legend(order(2 "IQR df=`df_str', {it:{&delta}}{sup:*}=`thres_str' (`fthres_str' Fisher transformed)") ///
                            pos(6) nobox region(color(none))) ///
                    xsize(16) ysize(9)
        graph export "`filesuffix'_hist.png", as(png) replace
    }
end   
* Function to save TMO results to dta file
cap program drop tmo_save
program define tmo_save,
    syntax, FILEsuffix(str)

    preserve
        clear

        qui set obs 1

        qui gen beta = beta
        qui gen orig_se = se
        qui gen tmo_se = tmo_se
        qui gen lb = lb
        qui gen ub = ub
        qui gen threshold = thres
        qui gen pct_ge_thres = offdP*100
        qui gen pct_ge_thres_nocl = offdPnocl*100
        qui gen T = T
        qui gen N_loc = N
        qui gen N = N_obs
        qui gen N_clust = N_clust
        qui gen N_outcomes = D
        qui gen dof = df
        qui gen finite_sample_dof = dof_adj
        qui gen df_r = df_r
        qui gen scpc_cv = ${scpc_cv}

        qui compress
        save "`filesuffix'_est.dta", replace
    restore

end