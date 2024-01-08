


*################################################################################################################################*
																* WALKING PACE *
*################################################################################################################################*

**# Bookmark #Analysis

**************************************
**Preparation*************************
**************************************
use "dataset_HD", replace
distinct id
drop if bcancer == 1
drop if bcvd == 1
tab cause_death
tab nmiss
drop fitness maxwl maxhr targethr ntrend durfit bcancer bcvd bckd bdiab n_22190_0_0 n_22191_0_0 n_22192_0_0 n_30000_0_0 mg100 nmiss probacc overacc wearacc

tab centre, missing
gen country = ""
replace country = "SCOT" if (centre == "Edinburgh" | centre == "Glasgow") & country == ""
replace country = "WAL"  if (centre == "Cardiff" | centre == "Swansea" | centre == "Wrexham") & country == ""
replace country = "ENG"  if (country == "") & country == ""
tab country, missing
order centre, before(country)

gen flag = 1 if dodn>date("31/7/2021", "DMY") & country == "SCOT"
replace dodn = date("31/7/2021", "DMY") if flag == 1
replace died = 0 if flag == 1
replace cause_death = "" if flag == 1
replace icd10 = "" if flag == 1
drop flag

gen flag = 1 if dodn>date("31/3/2016", "DMY") & country == "WAL"
replace dodn = date("31/3/2016", "DMY") if flag == 1
replace died = 0 if flag == 1
replace cause_death = "" if flag == 1
replace icd10 = "" if flag == 1													/*Unknown cause code with no icd10*/
drop flag

gen thosp = (t_cacvd-date0n)/365.24
label variable thosp "Time_to_hosp"
count if dodn == t_cacvd & (died == 1 & cacvd == 1)
list id dodn died cacvd t_cacvd if dodn == t_cacvd & (died == 1 & cacvd == 1) 	/*death same day of hospitalisation*/
tab cacvd, m
rename cacvd hosp
drop t_cacvd

gen tdied = (dodn-date0n)/365.24
label variable tdied "Time_to_death"
order died, before(tdied)
drop centre country date0 date0n dodn cause_death icd10 data_clean

sum tdied, d
tab hosp died, m
tab hosp died if tdied == thosp, m
tab wp, gen(dwp)

mdesc
describe
sort id
save "db_ms0", replace


**************************************
**Descriptive*************************
**************************************
clear all
use "db_ms0", clear
label define m_0009 0 "Women", modify
label define m_0009 1 "Men", modify

/*baseline tables by sex*/
baselinetable   											/*
*/	age0(cts tab("p50 (p25, p75)")) 						/*
*/  msbp(cts tab("p50 (p25, p75)") medianformat(%5.0f))		/*
*/  ldl(cts  tab("p50 (p25, p75)"))							/*
*/	bmi(cts tab("p50 (p25, p75)"))							/*
*/	tws(cts tab("p50 (p25, p75)")) 							/*
*/	tple(cts tab("p50 (p25, p75)")) 						/*
*/	smok(cat countformat(%15.0fc))							/*
*/	wp(cat countformat(%15.0fc))							/*
*/	, by(sex, total) countformat(%15.0fc) notable meanformat(%5.2f) sdformat(%5.1f) /*
*/	exportexcel("Results\Table_1", replace)

/*baseline tables by wp_sex*/
gen wp_sex = "wp" + string(wp) + "_sex" + string(sex)
baselinetable   											/*
*/	age0(cts tab("p50 (p25, p75)")) 						/*
*/  msbp(cts tab("p50 (p25, p75)") medianformat(%5.0f))		/*
*/  ldl(cts  tab("p50 (p25, p75)"))							/*
*/	bmi(cts tab("p50 (p25, p75)"))							/*
*/	tws(cts tab("p50 (p25, p75)")) 							/*
*/	tple(cts tab("p50 (p25, p75)")) 						/*
*/	smok(cat countformat(%15.0fc))							/*
*/	, by(wp_sex, total) countformat(%15.0fc) notable meanformat(%5.2f) sdformat(%5.2f) /*
*/	exportexcel("Results\Table_S1", replace)
drop wp_sex


*******************************
****Multistate model set-up****
*******************************
use "db_ms0", clear
gen tdn		= ""
replace tdn = "tdied>thosp" if tdied>thosp
replace tdn = "tdied=thosp" if tdied==thosp
replace tdn = "tdied<thosp" if tdied<thosp
tab tdn, m
sort tdn hosp died
groups tdn hosp died, sepby(tdn)

list id hosp thosp died tdied tdn 	 if hosp == 1 & died == 0 & tdn == "tdied=thosp"	
replace thosp = thosp - (0.5/365.24) if hosp == 1 & died == 0 & tdn == "tdied=thosp"	
list id hosp thosp died tdied tdn 	 if hosp == 1 & died == 1 & tdn == "tdied=thosp"
replace thosp = thosp - (0.5/365.24) if hosp == 1 & died == 1 & tdn == "tdied=thosp"	
list id hosp thosp died tdied if thosp == 0 | tdied == 0								
replace thosp = thosp + (0.5/365.24) if hosp == 1 & thosp == 0

replace tdn = "tdied>thosp" if tdied>thosp
replace tdn = "tdied=thosp" if tdied==thosp
replace tdn = "tdied<thosp" if tdied<thosp
tab tdn, m
sort tdn hosp died
groups tdn hosp died, sepby(tdn)
drop tdn

sort id
save "db_ms1", replace


***---------------------------------------------------------------------------------------------------------------------
***Rates - PROB, LOS, RMST [UNADJUSTED]---------------------------------------------------------------------------------
***---------------------------------------------------------------------------------------------------------------------
clear all
use "db_ms1", replace
sum tdied, d

preserve
keep if sex == 0
msset, id(id) states(hosp died) times(thosp tdied)
mat tmat = r(transmatrix)
mat list tmat
mat fmat = r(freqmatrix)
mat list fmat
matrix colnames fmat = "Healthy" "CVD/Cancer" "Died"
matrix rownames fmat = "Healthy" "CVD/Cancer" "Died"
mat list fmat
msboxes, transmat(tmat) id(id) xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) boxheight(0.3) boxwidth(0.3) ///
		 statenames("Women - Healthy [S1]" "CVD/Cancer [S2]" "Died [S3]") transnames("T1 (S1{&rarr}S2)" "T2 (S1{&rarr}S3)" "T3 (S2{&rarr}S3)")
restore

preserve
keep if sex == 1
msset, id(id) states(hosp died) times(thosp tdied)
mat tmat = r(transmatrix)
mat list tmat
mat fmat = r(freqmatrix)
mat list fmat
matrix colnames fmat = "Healthy" "CVD/Cancer" "Died"
matrix rownames fmat = "Healthy" "CVD/Cancer" "Died"
mat list fmat
msboxes, transmat(tmat) id(id) xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) boxheight(0.3) boxwidth(0.3) ///
		 statenames("Men - Healthy [S1]" "CVD/Cancer [S2]" "Died [S3]") transnames("T1 (S1{&rarr}S2)" "T2 (S1{&rarr}S3)" "T3 (S2{&rarr}S3)")
restore
graph close _all

msset, id(id) states(hosp died) times(thosp tdied)
mat tmat = r(transmatrix)
mat list tmat
mat fmat = r(freqmatrix)
mat list fmat
matrix colnames fmat = "Healthy" "CVD/Cancer" "Died"
matrix rownames fmat = "Healthy" "CVD/Cancer" "Died"
mat list fmat
replace age0 = age0 + _start if _from == 2
stset _stop, enter(_start) failure(_status==1)

preserve
clear
tempfile results
save `results', emptyok replace
restore

forvalues s = 0/1 {
	
	mata: mata clear
	
	preserve
	
	qui keep if sex == `s'
	
	/*First transition*/
	qui strate wp if _trans1 == 1, per(100) output("Results\RateSex`s'_T1", replace)
	qui stmerlin dwp2 dwp3 age0 if _trans1 == 1, dist(rp) df(4) nolog	
	qui estimate store t1

	/*Second transition*/
	qui strate wp if _trans2 == 1, per(100) output("Results\RateSex`s'_T2", replace)
	qui stmerlin dwp2 dwp3 age0 if _trans2 == 1, dist(rp) df(4) nolog	
	qui estimate store t2

	/*Third transition*/
	qui strate wp if _trans3 == 1, per(100) output("Results\RateSex`s'_T3", replace)
	qui stmerlin dwp2 dwp3 age0 if _trans3 == 1, dist(rp) df(4) nolog	
	qui estimate store t3
	
	/*Numerical integration*/
	qui tempfile mainfile
	qui save `mainfile', replace
	
	forvalues a = 55(10)75 {															
		qui gen tvar = 10 in 1														
		qui predictms, transmatrix(tmat) models(t1 t2 t3) timevar(tvar) ///
				   at1(dwp2 0 dwp3 0 age0 `a') 			///
				   at2(dwp2 1 dwp3 0 age0 `a') 			///
				   at3(dwp2 0 dwp3 1 age0 `a') 			///
				   probability rmst los difference ci
		qui gen agep = `a'
		qui gen sexp = `s'
		qui keep sexp agep tvar _prob* _rmst* _los* _diff* 
		qui drop if tvar ==.
		qui append using `results'
		qui save `results', replace
		qui use `mainfile', clear
		di " -- Sex = `s' | Age = `a' -- $S_TIME $S_DATE -- "
	}
	restore
}
use `results', clear
gen model = "Unadjusted"
order model sexp agep tvar, first
save "Results\Unadjusted", replace

clear
tempfile rates
save `rates', emptyok replace
forvalues s = 0/1 {
	forvalues t = 1/3 {
		qui use  "Results\RateSex`s'_T`t'", clear
		qui gen Sex = `s'
		qui gen Transition = "T`t'"
		qui append using `rates'
		qui save `rates', replace
	}
}
use `rates', clear
gen 	States = "Healthy [S1] to CVD/Cancer [S2]" 	if Transition == "T1"
replace States = "Healthy [S1] to Died [S3]" 	 	if Transition == "T2"
replace States = "CVD/Cancer [S2] to Died [S3]" 	if Transition == "T3"
renames _D _Y _Rate _Lower _Upper \ Events PYRs Rate100 LB95 UB95
tostring Sex, replace
replace Sex = "Women" if Sex == "0"
replace Sex = "Men" if Sex == "1"
rename wp WalkingPace
compress
order Sex WalkingPace Transition States, first
export delimited using "Results\Table_S2.csv", datafmt replace

forvalues s = 0/1 {
	forvalues t = 1/3 {
		qui erase  "Results\RateSex`s'_T`t'.dta"
	}
}


***---------------------------------------------------------------------------------------------------------------------
***PROB, LOS, RMST [ADJUSTED/CONDITIONAL]-------------------------------------------------------------------------------
***---------------------------------------------------------------------------------------------------------------------
use "db_ms1", replace
tab smok, gen(dsmok)

msset, id(id) states(hosp died) times(thosp tdied)
mat tmat = r(transmatrix)
replace age0 = age0 + _start if _from == 2
stset _stop, enter(_start) failure(_status==1)

preserve
clear
tempfile results
save `results', emptyok replace
restore

forvalues s = 0/1 {
	
	mata: mata clear
	
	preserve
	
	qui keep if sex == `s'
	
	foreach var of varlist msbp ldl tws tple bmi {
	qui sum `var', meanonly
	qui local `var'_m = r(mean)
	}
	
	/*First transition*/
	qui stmerlin dwp2 dwp3 age0 msbp ldl tws tple dsmok2 dsmok3 bmi if _trans1 == 1, dist(rp) df(4) nolog	
	qui estimate store t1

	/*Second transition*/
	qui stmerlin dwp2 dwp3 age0 msbp ldl tws tple dsmok2 dsmok3 bmi if _trans2 == 1, dist(rp) df(4) nolog	
	qui estimate store t2

	/*Third transition*/
	qui stmerlin dwp2 dwp3 age0 msbp ldl tws tple dsmok2 dsmok3 bmi if _trans3 == 1, dist(rp) df(4) nolog	
	qui estimate store t3
	
	/*Numerical integration*/
	qui tempfile mainfile
	qui save `mainfile', replace
	
	forvalues a = 55(10)75 {																								
		qui gen tvar = 10 in 1																							
		qui predictms, transmatrix(tmat) models(t1 t2 t3) timevar(tvar) ///
				   at1(dwp2 0 dwp3 0 age0 `a' msbp `msbp_m' ldl `ldl_m' tws `tws_m' tple `tple_m' dsmok2 0 dsmok3 0 bmi `bmi_m') 	/// /*for most common, no smoker*/
				   at2(dwp2 1 dwp3 0 age0 `a' msbp `msbp_m' ldl `ldl_m' tws `tws_m' tple `tple_m' dsmok2 0 dsmok3 0 bmi `bmi_m') 	///
				   at3(dwp2 0 dwp3 1 age0 `a' msbp `msbp_m' ldl `ldl_m' tws `tws_m' tple `tple_m' dsmok2 0 dsmok3 0 bmi `bmi_m')	///
				   probability rmst los difference ci
		qui gen agep  = `a'
		qui gen sexp  = `s'
		qui gen msbpp = `msbp_m'
		qui gen ldlp  = `ldl_m'
		qui gen twsp  = `tws_m'
		qui gen tplep = `tple_m'
		qui gen smokp = "Never"
		qui gen bmip  = `bmi_m'
		qui keep agep-bmip tvar _prob* _rmst* _los* _diff* 
		qui drop if tvar ==.
		qui append using `results'
		qui save `results', replace
		qui use `mainfile', clear
		di " -- Sex = `s' | Age = `a' -- $S_TIME $S_DATE -- "
	}
	restore
}
use `results', clear
gen model = "Adjusted"
order model sexp agep msbpp-bmip tvar, first
save "Results\Adjusted", replace



**# Bookmark #Graphs

************************
****GRAPHS**************
************************

set scheme white_tableau

****ABSOLUTE VALUES****
use "Unadjusted", clear
append using "Adjusted"
order model sexp agep msbpp-bmip tvar, first
drop msbpp-bmip _diff*

*Probabilitites and LOS - preparation*
preserve
drop _rmst*
reshape long _prob_ _los_, i(model sexp agep tvar) j(metric) string
gen metric1 = ""
replace metric1 = subinstr(metric, "at1_1_1", "[S1] Healthy [WP1] Slow", .) 		if substr(metric, 1, 7) == "at1_1_1"
replace metric1 = subinstr(metric, "at1_1_2", "[S2] CVD/Cancer [WP1] Slow", .) 		if substr(metric, 1, 7) == "at1_1_2"
replace metric1 = subinstr(metric, "at1_1_3", "[S3] Died [WP1] Slow" , .) 			if substr(metric, 1, 7) == "at1_1_3"
replace metric1 = subinstr(metric, "at2_1_1", "[S1] Healthy [WP2] Average", .) 		if substr(metric, 1, 7) == "at2_1_1"
replace metric1 = subinstr(metric, "at2_1_2", "[S2] CVD/Cancer [WP2] Average", .)	if substr(metric, 1, 7) == "at2_1_2"
replace metric1 = subinstr(metric, "at2_1_3", "[S3] Died [WP2] Average" , .) 		if substr(metric, 1, 7) == "at2_1_3"
replace metric1 = subinstr(metric, "at3_1_1", "[S1] Healthy [WP3] Brisk", .) 		if substr(metric, 1, 7) == "at3_1_1"
replace metric1 = subinstr(metric, "at3_1_2", "[S2] CVD/Cancer [WP3] Brisk", .) 	if substr(metric, 1, 7) == "at3_1_2"
replace metric1 = subinstr(metric, "at3_1_3", "[S3] Died [WP3] Brisk" , .) 			if substr(metric, 1, 7) == "at3_1_3"
drop metric
split metric1, p("")
split metric14, p("_")
replace metric142 = "_e" if metric142 == ""
replace metric142 = "_" + metric142 if metric142 != "_e"
replace metric1 = metric11 + " " + metric12 + " " + metric13 + " " + metric141
drop metric14
renames metric1 metric11 metric12 metric13 metric141 metric142 _prob_ _los_ \ group staten state pacen pace estimate PROB LOS
reshape wide PROB LOS, i(model sexp agep tvar group) j(estimate) string
sort model sexp agep pacen staten tvar
order PROB* LOS*, after(group)
renames PROB_e-LOS_uci \ estPROB lciPROB uciPROB estLOS lciLOS uciLOS
reshape long est lci uci, i(model sexp agep tvar group) j(metric) string
tostring sexp, replace
replace sexp = "Women" if sexp == "0"
replace sexp = "Men" if sexp == "1"
renames sexp agep tvar \ sex age time
sencode sex, gsort(-sex) replace
sencode pace, gsort(pacen) replace
sort metric model sex age pacen staten
compress
tempfile abs_prob_los
save `abs_prob_los', replace
restore

*Graph LOS*
use `abs_prob_los', replace
keep if metric == "LOS"
drop time lci uci
sencode state, gsort(staten) replace
drop if staten == "[S3]"
tempfile los
save `los', replace

graph bar (asis) est if model == "Adjusted", over(state, gap(200)) over(pace) asyvars stack ///
	  ylabel(, gmax) by(sex age, cols(3) note("")) subtitle(, size(small))					///
	  legend(rows(1) size(small)) ytitle("Length of stay, yrs", size(medsmall)) 			///
	  xsize(6) ysize(4) name("Fig1", replace) nodraw

graph bar (asis) est if model == "Unadjusted", over(state, gap(200)) over(pace) asyvars stack 	///
	  ylabel(, gmax) by(sex age, cols(3) note("")) subtitle(, size(small))						///
	  legend(rows(1) size(small)) ytitle("Length of stay, yrs", size(medsmall)) 				///
	  xsize(6) ysize(4) name("FigS3", replace) nodraw


****DIFFERENCES VALUES****	  
use "Unadjusted", clear
append using "Adjusted"
order model sexp agep msbpp-bmip tvar, first
drop msbpp-bmip _prob* _los* _rmst* _diff_prob*

*LOS - preparation*
preserve
drop _diff_rmst*
reshape long _diff_los_, i(model sexp agep tvar) j(metric) string
gen metric1 = ""
replace metric1 = subinstr(metric, "at2_1_1", "[S1] Healthy [WP2 vs WP1] Average vs Slow", .) 		if substr(metric, 1, 7) == "at2_1_1"
replace metric1 = subinstr(metric, "at2_1_2", "[S2] CVD/Cancer [WP2 vs WP1] Average vs Slow", .) 	if substr(metric, 1, 7) == "at2_1_2"
replace metric1 = subinstr(metric, "at2_1_3", "[S3] Died [WP2 vs WP1] Average vs Slow" , .) 		if substr(metric, 1, 7) == "at2_1_3"
replace metric1 = subinstr(metric, "at3_1_1", "[S1] Healthy [WP3 vs WP1] Brisk vs Slow", .) 		if substr(metric, 1, 7) == "at3_1_1"
replace metric1 = subinstr(metric, "at3_1_2", "[S2] CVD/Cancer [WP3 vs WP1] Brisk vs Slow", .) 		if substr(metric, 1, 7) == "at3_1_2"
replace metric1 = subinstr(metric, "at3_1_3", "[S3] Died [WP3 vs WP1] Brisk vs Slow" , .) 			if substr(metric, 1, 7) == "at3_1_3"
drop metric
split metric1, p("_")
replace metric1 = metric11
drop metric11
replace metric12 = "_e" if metric12 == ""
replace metric12 = "_" + metric12 if metric12 != "_e"
rename metric1 metric
split metric, p("")
gen dwpn = metric3 + " " + metric4 + " " + metric5
gen dwp = metric6 + " " + metric7 + " " + metric8
drop metric3-metric8 tvar
renames _diff_los_ metric metric12 metric1 metric2 \ DLOS comparison estimate staten state
reshape wide DLOS, i(model sexp agep comparison) j(estimate) string
renames DLOS_e DLOS_lci DLOS_uci \ est lci uci
gen metric = "DLOS", before(est)
sort model sexp agep dwpn staten
tostring sexp, replace
replace sexp = "Women" if sexp == "0"
replace sexp = "Men" if sexp == "1"
renames sexp agep \ sex age 
compress
tempfile diff_los
save `diff_los', replace
restore

*RMST - preparation*
drop _diff_los* tvar
reshape long _diff_rmst_, i(model sexp agep) j(metric) string
replace metric = subinstr(metric, "_1", "", .)									
replace metric = subinstr(metric, "at2", "[WP2 vs WP1] Average vs Slow", .) 	
replace metric = subinstr(metric, "at3", "[WP3 vs WP1] Brisk vs Slow", .)
renames metric _diff_rmst_ \ comparison DRMST
split comparison, p("_")
replace comparison = comparison1
replace comparison2 = "_e" if comparison2 == ""
replace comparison2 = "_" + comparison2 if comparison2 != "_e"
split comparison1, p("] ")
replace comparison11 = comparison11 + "]"
drop comparison1
renames comparison11 comparison12 \ dwpn dwp
reshape wide DRMST, i(model sexp agep comparison) j(comparison2) string
renames DRMST_e DRMST_lci DRMST_uci \ est lci uci
gen metric = "DRMST", before(est)
gen staten = "[S1|S2]", before(dwpn)
gen state = "NotDied", before(dwpn)
sort model sexp agep dwpn
tostring sexp, replace
replace sexp = "Women" if sexp == "0"
replace sexp = "Men" if sexp == "1"
renames sexp agep \ sex age 
compress
tempfile diff_rmst
save `diff_rmst', replace

*Graphs LOS & RMST*
use `diff_los', replace
append using `diff_rmst'
order metric, first
sencode state, replace
sencode sex, gsort(-sex) replace
foreach var of varlist est-uci {
	replace `var' = `var'*12
}
drop if staten == "[S3]"
sort model sex age state dwpn
tempfile diffs
save `diffs', replace

egen float grp = group(dwp)
sdecode sex, replace

foreach s in Women Men {
	foreach a in Adjusted Unadjusted {
		forestplot est lci uci if model == "`a'" & sex == "`s'", lcols(age dwp state) nonames nooverall nosubgroup dp(1) 	///
		xlabel(-12 -8 -4 0 4 8 12 16 18 20, force) classic boxscale(85) astext(50) ciopts(lwidth(vthin)) plotid(state)		///
		box1opts(mcolor(blue)) ci1opts(lcolor(blue)) box2opts(mcolor(orange)) ci2opts(lcolor(orange)) xtitle("Months")		///
		leftjustify name("`s'_`a'", replace) title("`s'", size(medsmall)) nodraw
	}
}

graph combine Women_Adjusted   Men_Adjusted,   cols(2) nocopies name("Fig2", replace) 
graph combine Women_Unadjusted Men_Unadjusted, cols(2) nocopies name("FigS4", replace) 
graph close _all


**********************************************
****ALL GRAPHS********************************
**********************************************

*[Figure 1 - Adjusted LOS]*
use `los', replace

*[Figure 2 - Adjusted difference in LOS and RMST]*
use `diffs', replace

*[Figure S1 - Flowchart of study participants]*
*[Figure S2 - Transition diagram]*
*[Figure S3 - Unadjusted LOS]*
use `los', replace

*[Figure S4 - Unadjusted difference in LOS and RMST]*
use `diffs'

foreach g in Fig1 Fig2 FigS3 FigS4 {
	graph save "`g'" "`g'.gph", replace
}
graph close _all



*################################################################################################################################*
													* WALKING PACE --- USING FIRST CAUSE OF HOSPITALISATION*
*################################################################################################################################*

**# Bookmark #ONLY FIRST ICD
use "dataset_HD", replace
keep id date0n
tempfile dtn
save `dtn', replace

use "db_ms0", clear
merge 1:1 id using "`dtn'", update
keep if _merge == 3
drop _merge hosp thosp

merge 1:1 id using "db_hospitalisation_icd1", update
keep if _merge == 3
drop _merge

gen thosp_icd1 = (t_cacvd_icd1-date0n)/365.24
drop date0n t_cacvd_icd1 data_clean
rename cacvd_icd1 hosp_icd1

mdesc
describe
save "db_ms0_icd1", replace


********************************************************
**Set-up transitions [see above same part for details]**
********************************************************
use "db_ms0_icd1", clear

gen tdn		= ""
replace tdn = "tdied>thosp" if tdied>thosp_icd1
replace tdn = "tdied=thosp" if tdied==thosp_icd1
replace tdn = "tdied<thosp" if tdied<thosp_icd1
tab tdn, m
sort tdn hosp died
groups tdn hosp_icd1 died, sepby(tdn)

replace thosp_icd1 = thosp_icd1 - (0.5/365.24) if hosp_icd1 == 1 & died == 0 & tdn == "tdied=thosp"	
replace thosp_icd1 = thosp_icd1 - (0.5/365.24) if hosp_icd1 == 1 & died == 1 & tdn == "tdied=thosp"
replace thosp_icd1 = thosp_icd1 + (0.5/365.24) if hosp_icd1 == 1 & thosp_icd1 == 0					

replace tdn = "tdied>thosp" if tdied>thosp_icd1
replace tdn = "tdied=thosp" if tdied==thosp_icd1
replace tdn = "tdied<thosp" if tdied<thosp_icd1
tab tdn, m
sort tdn hosp died
groups tdn hosp died, sepby(tdn)
drop tdn

order died tdied hosp_icd1 thosp_icd1, last
sort id
save "db_ms1_icd1", replace


************************
****Multistate model****
************************
use "db_ms1_icd1", replace

preserve
keep if sex == 0
msset, id(id) states(hosp died) times(thosp tdied)
mat tmat = r(transmatrix)
mat list tmat
mat fmat = r(freqmatrix)
mat list fmat
matrix colnames fmat = "Healthy" "CVD/Cancer" "Died"
matrix rownames fmat = "Healthy" "CVD/Cancer" "Died"
mat list fmat
msboxes, transmat(tmat) id(id) xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) boxheight(0.3) boxwidth(0.3) ///
		 statenames("Women - Healthy [S1]" "CVD/Cancer [S2]" "Died [S3]") transnames("T1 (S1{&rarr}S2)" "T2 (S1{&rarr}S3)" "T3 (S2{&rarr}S3)")
restore

preserve
keep if sex == 1
msset, id(id) states(hosp died) times(thosp tdied)
mat tmat = r(transmatrix)
mat list tmat
mat fmat = r(freqmatrix)
mat list fmat
matrix colnames fmat = "Healthy" "CVD/Cancer" "Died"
matrix rownames fmat = "Healthy" "CVD/Cancer" "Died"
mat list fmat
msboxes, transmat(tmat) id(id) xvalues(0.2 0.7 0.45) yvalues(0.7 0.7 0.2) boxheight(0.3) boxwidth(0.3) ///
		 statenames("Men - Healthy [S1]" "CVD/Cancer [S2]" "Died [S3]") transnames("T1 (S1{&rarr}S2)" "T2 (S1{&rarr}S3)" "T3 (S2{&rarr}S3)")
restore
graph close _all


***---------------------------------------------------------------------------------------------------------------------
***Probabiltites, LOS, RMST [ADJUSTED/CONDITIONAL]----------------------------------------------------------------------
***---------------------------------------------------------------------------------------------------------------------
use "db_ms1_icd1", replace
tab smok, gen(dsmok)
msset, id(id) states(hosp_icd1 died) times(thosp_icd1 tdied)
mat tmat = r(transmatrix)
replace age0 = age0 + _start if _from == 2
stset _stop, enter(_start) failure(_status==1)

preserve
clear
tempfile results
save `results', emptyok replace
restore

forvalues s = 0/1 {
	
	mata: mata clear
	
	preserve
	
	qui keep if sex == `s'
	
	foreach var of varlist msbp ldl tws tple bmi {
	qui sum `var', meanonly
	qui local `var'_m = r(mean)
	}
	
	/*First transition*/
	qui stmerlin dwp2 dwp3 age0 msbp ldl tws tple dsmok2 dsmok3 bmi if _trans1 == 1, dist(rp) df(4) nolog	
	qui estimate store t1

	/*Second transition*/
	qui stmerlin dwp2 dwp3 age0 msbp ldl tws tple dsmok2 dsmok3 bmi if _trans2 == 1, dist(rp) df(4) nolog	
	qui estimate store t2

	/*Third transition*/
	qui stmerlin dwp2 dwp3 age0 msbp ldl tws tple dsmok2 dsmok3 bmi if _trans3 == 1, dist(rp) df(4) nolog	
	qui estimate store t3
	
	/*Numerical integration*/
	qui tempfile mainfile
	qui save `mainfile', replace
	
	forvalues a = 55(10)75 {																								
		qui gen tvar = 10 in 1																							
		qui predictms, transmatrix(tmat) models(t1 t2 t3) timevar(tvar) ///
				   at1(dwp2 0 dwp3 0 age0 `a' msbp `msbp_m' ldl `ldl_m' tws `tws_m' tple `tple_m' dsmok2 0 dsmok3 0 bmi `bmi_m') 	/// /*for most common, no smoker*/
				   at2(dwp2 1 dwp3 0 age0 `a' msbp `msbp_m' ldl `ldl_m' tws `tws_m' tple `tple_m' dsmok2 0 dsmok3 0 bmi `bmi_m') 	///
				   at3(dwp2 0 dwp3 1 age0 `a' msbp `msbp_m' ldl `ldl_m' tws `tws_m' tple `tple_m' dsmok2 0 dsmok3 0 bmi `bmi_m')	///
				   probability rmst los difference ci
		qui gen agep  = `a'
		qui gen sexp  = `s'
		qui gen msbpp = `msbp_m'
		qui gen ldlp  = `ldl_m'
		qui gen twsp  = `tws_m'
		qui gen tplep = `tple_m'
		qui gen smokp = "Never"
		qui gen bmip  = `bmi_m'
		qui keep agep-bmip tvar _prob* _rmst* _los* _diff* 
		qui drop if tvar ==.
		qui append using `results'
		qui save `results', replace
		qui use `mainfile', clear
		di " -- Sex = `s' | Age = `a' -- $S_TIME $S_DATE -- "
	}
	restore
}
use `results', clear
gen model = "Adjusted_ICD1"
order model sexp agep msbpp-bmip tvar, first
save "Results\Adjusted_ICD1", replace


**********************************************
********************GRAPHS********************
**********************************************
use "Adjusted_ICD1", clear
order model sexp agep msbpp-bmip tvar, first
drop msbpp-bmip _prob* _los* _rmst* _diff_prob*

*LOS - preparation*
preserve
drop _diff_rmst*
reshape long _diff_los_, i(model sexp agep tvar) j(metric) string
gen metric1 = ""
replace metric1 = subinstr(metric, "at2_1_1", "[S1] Healthy [WP2 vs WP1] Average vs Slow", .) 		if substr(metric, 1, 7) == "at2_1_1"
replace metric1 = subinstr(metric, "at2_1_2", "[S2] CVD/Cancer [WP2 vs WP1] Average vs Slow", .) 	if substr(metric, 1, 7) == "at2_1_2"
replace metric1 = subinstr(metric, "at2_1_3", "[S3] Died [WP2 vs WP1] Average vs Slow" , .) 		if substr(metric, 1, 7) == "at2_1_3"
replace metric1 = subinstr(metric, "at3_1_1", "[S1] Healthy [WP3 vs WP1] Brisk vs Slow", .) 		if substr(metric, 1, 7) == "at3_1_1"
replace metric1 = subinstr(metric, "at3_1_2", "[S2] CVD/Cancer [WP3 vs WP1] Brisk vs Slow", .) 		if substr(metric, 1, 7) == "at3_1_2"
replace metric1 = subinstr(metric, "at3_1_3", "[S3] Died [WP3 vs WP1] Brisk vs Slow" , .) 			if substr(metric, 1, 7) == "at3_1_3"
drop metric
split metric1, p("_")
replace metric1 = metric11
drop metric11
replace metric12 = "_e" if metric12 == ""
replace metric12 = "_" + metric12 if metric12 != "_e"
rename metric1 metric
split metric, p("")
gen dwpn = metric3 + " " + metric4 + " " + metric5
gen dwp = metric6 + " " + metric7 + " " + metric8
drop metric3-metric8 tvar
renames _diff_los_ metric metric12 metric1 metric2 \ DLOS comparison estimate staten state
reshape wide DLOS, i(model sexp agep comparison) j(estimate) string
renames DLOS_e DLOS_lci DLOS_uci \ est lci uci
gen metric = "DLOS", before(est)
sort model sexp agep dwpn staten
tostring sexp, replace
replace sexp = "Women" if sexp == "0"
replace sexp = "Men" if sexp == "1"
renames sexp agep \ sex age 
compress
tempfile diff_los
save `diff_los', replace
restore

*RMST - preparation*
drop _diff_los* tvar
reshape long _diff_rmst_, i(model sexp agep) j(metric) string
replace metric = subinstr(metric, "_1", "", .)									
replace metric = subinstr(metric, "at2", "[WP2 vs WP1] Average vs Slow", .) 	
replace metric = subinstr(metric, "at3", "[WP3 vs WP1] Brisk vs Slow", .)
renames metric _diff_rmst_ \ comparison DRMST
split comparison, p("_")
replace comparison = comparison1
replace comparison2 = "_e" if comparison2 == ""
replace comparison2 = "_" + comparison2 if comparison2 != "_e"
split comparison1, p("] ")
replace comparison11 = comparison11 + "]"
drop comparison1
renames comparison11 comparison12 \ dwpn dwp
reshape wide DRMST, i(model sexp agep comparison) j(comparison2) string
renames DRMST_e DRMST_lci DRMST_uci \ est lci uci
gen metric = "DRMST", before(est)
gen staten = "[S1|S2]", before(dwpn)
gen state = "NotDied", before(dwpn)
sort model sexp agep dwpn
tostring sexp, replace
replace sexp = "Women" if sexp == "0"
replace sexp = "Men" if sexp == "1"
renames sexp agep \ sex age 
compress
tempfile diff_rmst
save `diff_rmst', replace

*Graphs LOS & RMST*
use `diff_los', replace
append using `diff_rmst'
order metric, first
sencode state, replace
sencode sex, gsort(-sex) replace
foreach var of varlist est-uci {
	replace `var' = `var'*12
}
drop if staten == "[S3]"
sort model sex age state dwpn
tempfile diffs_icd1
save `diffs_icd1', replace

egen float grp = group(dwp)
sdecode sex, replace

foreach s in Women Men {
	forestplot est lci uci if sex == "`s'", lcols(age dwp state) nonames nooverall nosubgroup dp(1) 	            ///
	xlabel(-6 -3 0 3 6 9 12, force) classic boxscale(85) astext(50) ciopts(lwidth(vthin)) plotid(state)	            ///
	box1opts(mcolor(blue)) ci1opts(lcolor(blue)) box2opts(mcolor(orange)) ci2opts(lcolor(orange)) xtitle("Months")	///
	leftjustify name("`s'_icd1", replace) title("`s'_icd1", size(medsmall)) nodraw
}

graph combine Women_icd1 Men_icd1, cols(2) nocopies name("FigS6", replace) nodraw


**********************************************
****ALL GRAPHS********************************
**********************************************

*[Figure S5 - Transition ICD first only]*

*[Figure S6 - Adjusted difference in LOS and RMST using ICD fist only]*
use `diffs_icd1', replace

graph save "FigS6" "FigS6.gph", replace
graph close _all