* ---------------------------------------------------------------------------- *
* --------			 		FCI CONSTRUCTION FOR CHILE 				---------- *
* --------					Coded by: Hriday Karnani 				---------- *
* Date: October 2022
* ---------------------------------------------------------------------------- *

* Introduction to the code:
* In this code I create a Financial Conditions Index (FCI) based on a set of 
* financial variables. The methodology is a simple PCA, but using a 
* de-factored version of these financial variables. Later, I compare and
* estimate the impact of the FCI on economic activity with simple OLS and 
* quantile regressions. Finally, I estimate the out of sample performance of
* this FCI.

clear
set more off

* Define working directory
global maindir "C:/Users/Hriday/Dropbox/FCI CL/"
global construction "${maindir}FCI_Construction/"
global figs "${construction}Figs/"
global tabs "${construction}Tabs/"
global data"${construction}Data/"

cd "${maindir}"

* This dofile uses grqreg, corrtex, frmttable and st0207 packages
* ssc install grqreg
* ssc install corrtex
* ssc install frmttable
* ssc install st0207

* Open dataset
import excel "${maindir}FCI Series.xlsx", sheet("Agg_Data") firstrow clear

* Gen Nominal Exchange Rate as a m-o-m variation
gen TCNmom = TCN[_n]/TCN[_n-1] -1
* Gen Imacec growth (economic activity), m-o-m
gen Imacec_g = log(Imacec_Des[_n]) - log(Imacec_Des[_n-1])
drop TCN
* Manipulate date format in order to make regression analysis easier later
gen Month = tm(1990m1) + _n-1
format %tm Month
drop Date 
rename Month Date 
order Date
tsset Date

* Simple sanity checks of the data: Correlation matrix
* Export
corrtex Corp_Sov_Spread MPR Term_Spread EMBI IPSA_Vol G_Credit_GDP Chicago_FCI VIX TCNmom, file(Corr_Matrix) replace 


* Defactor each financial variables using Macro variables as regressors:
* Standardize variables (although PCA funcion does it on its own)
foreach v in Corp_Sov_Spread MPR Term_Spread EMBI IPSA_Vol G_Credit_GDP Chicago_FCI VIX TCNmom {
	egen mean = mean(`v')
	egen sd = sd(`v')
	gen std_`v' = (`v'-mean)/sd
	drop mean sd 
	reg std_`v' Inflation L1.Inflation L2.Inflation GDP_Monthly L1.GDP_Monthly L2.GDP_Monthly, r
	predict resid_`v', residuals
	drop `v' std_`v'
	rename resid_`v' `v'
}


* PC Analysis
pca Corp_Sov_Spread MPR Term_Spread EMBI IPSA_Vol G_Credit_GDP Chicago_FCI VIX TCNmom, comp(3)
* Export Eigenvectors
matrix eigenvectors = e(L)
frmttable using "${tabs}Eigenvectors_Haztius.tex", statmat(eigenvectors) sfmt(g) sdec(3) ctitle("" "Eigenvectors 1st PC" "2PC" "3PC") fragment tex  hline(11000000001) vlines(01001) replace

* Create aux variables to estimate the contribution of each variable on the
* total FCI. Multiply eigenvector with the standardize value of the financial
* variable. This will give us the contribution of each variable on the FCI.
* in each iteration I sum the contribution of each variable in order to makes
* a "stacked" plot that nicely shows the contribution of each variable.

local i = 1
foreach v in Corp_Sov_Spread MPR Term_Spread EMBI IPSA_Vol G_Credit_GDP Chicago_FCI VIX TCNmom {
	gen pc1_`i' = eigenvectors[`i',1]*`v'	
	local j = `i'-1
	if `i'==1{
		gen aux_pc1_`i' = pc1_`i'
	}
	else {
		gen aux_pc1_`i' = pc1_`i' + aux_pc1_`j'
	}
	label variable aux_pc1_`i' `v'
	local i = `i'+1
}


* Check eigenvalues
screeplot
* Export
graph export "${figs}Eigenvalues.eps",replace
* Create index with PC1
predict pc1 pc2 pc3, score
label variable pc1 "FCI"


* Plot the series of the eigenvector gained from PCA, which can be interpreted as FCI
line pc1 Date,title("Financial Conditions Index 1PC") title("Financial Conditions Index for Chile") note("Note: A higher value of the index indicates tighter financial conditions") scheme(s1mono)
graph export "${figs}FCI_Haztius.eps",replace

* Plot the decomposition of the FCI: Watch the contribution of each variable 
* on different periods of time.
twoway area aux_pc1_9 aux_pc1_8 aux_pc1_7 aux_pc1_6 aux_pc1_5 aux_pc1_4 aux_pc1_3 aux_pc1_2 aux_pc1_1  Date, xla(360(60)750) scheme(s1color) title("FCI Decomposition - Chile") xtitle("Date") ytitle("FCI") note("Note: A higher value of the index indicates tighter financial conditions") legend( col(3)) 
graph export "${figs}FCI_Decomposition_Haztius.eps",replace

* Compare the FCI with economic activity on a plot
gen pc1_neg = -pc1
* I flipped the FCI in order to compare better Activity growth and our FCI
* FCI and Activity growth in the same plot
twoway (line Imacec_g Date, yaxis(1)) (line pc1_neg Date, yaxis(2)), scheme(s1color) legend(label(1 Imacec Growth (LHS)) label(2 FCI (RHS)))  title("FCI and Imacec Growth in Chile") ytitle("Imacec Growth") ytitle("FCI", axis(2)) note("Note: I flipped the FCI direction in order to make variables comparable.")
graph export "${figs}FCI_Imacecg_MOMDes_Haztius.eps",replace

* Label variables for regressions and gen lags and forwards for quantile reg
* (quan reg function doesn't allow L. or F. commands)
label variable Imacec_g "Imacec growth t"
gen l1_Imacec_g = Imacec_g[_n-1]
label variable l1_Imacec_g "Imacec growth t-1"
gen Imacec_g_f1 = Imacec_g[_n+1]
label variable Imacec_g_f1 "Imacec growth t+1"
gen Imacec_g_f2 = Imacec_g[_n+2]
label variable Imacec_g_f1 "Imacec growth t+2"
gen Imacec_g_f3 = Imacec_g[_n+3]
label variable Imacec_g_f1 "Imacec growth t+3"
gen Imacec_g_f5 = Imacec_g[_n+5]
label variable Imacec_g_f1 "Imacec growth t+5"
gen Imacec_g_f6 = Imacec_g[_n+6]
label variable Imacec_g_f1 "Imacec growth t+6"

* ---------------------------------------------------------------------------- *
* Preliminary regressions
* ---------------------------------------------------------------------------- *
* Compare FCI with leads and present economic activity
eststo M1: reg Imacec_g pc1,r
eststo M2: reg Imacec_g L.Imacec_g pc1,r
predict hat_imacec_g 
eststo M3: reg F1.Imacec_g L.Imacec_g pc1,r
eststo M4: reg F2.Imacec_g L.Imacec_g pc1,r
eststo M5: reg F3.Imacec_g L.Imacec_g pc1,r
eststo M6: reg F6.Imacec_g L.Imacec_g pc1,r
* Export
esttab M1 M2 M3 M4 M5 M6 using "${tabs}OLS_MOMDes_Haztius.tex", replace nogaps b(%-9.3f) se(%-9.3f)  starlevel("*" 0.10 "**" 0.05 "***" 0.01) stats(N r2, fmt(%12.0gc %04.2f)  label("Obs." "R-Squared")) label nonotes nonumbers style(tex) mtitles("ImacecG t" "ImacecG t" "ImacecG t+1" "ImacecG t+2" "ImacecG t+3" "ImacecG t+6") 
drop _est*

* Include lagged economic activity as a regressor (trying to solve some 
* endogeneity issues, obviously not completely)
eststo M1: reg Imacec_g pc1,r
eststo M2: reg Imacec_g L.Imacec_g pc1,r
eststo M3: reg F1.Imacec_g F1.L.Imacec_g pc1,r
eststo M4: reg F2.Imacec_g F2.L.Imacec_g pc1,r
eststo M5: reg F3.Imacec_g F3.L.Imacec_g pc1,r
eststo M6: reg F6.Imacec_g F6.L.Imacec_g pc1,r

esttab M1 M2 M3 M4 M5 M6 using "${tabs}OLS_MOMDest-1_Haztius.tex", replace b(%-9.3f) se(%-9.3f)  starlevel("*" 0.10 "**" 0.05 "***" 0.01) stats(N r2, fmt(%12.0gc %04.2f)  label("Obs." "R-Squared")) label nogaps nonotes nonumbers style(tex) mtitles("ImacecG t" "ImacecG t" "ImacecG t+1" "ImacecG t+2" "ImacecG t+3" "ImacecG t+6") 

drop _est*

* ---------------------------------------------------------------------------- *
* Quantile regressions
* -----------------------------------------------------------------------------*
* Quantile regressions for quantiles 5, 10, 50 and 90. Different versions with 
* leads of economic activity, with and without lagged economic activity as
* regressor
foreach quan of numlist 5 10 50 90 {
	eststo Q1M1: quietly qreg Imacec_g pc1, vce(robust) quantile(.`quan')
	eststo Q1M2: quietly qreg Imacec_g l1_Imacec_g pc1,vce(robust) quantile(.`quan')
	eststo Q1M3: quietly qreg Imacec_g_f1 l1_Imacec_g pc1,vce(robust) quantile(.`quan')
	eststo Q1M4: quietly qreg Imacec_g_f2 l1_Imacec_g pc1,vce(robust) quantile(.`quan')
	eststo Q1M5: quietly qreg Imacec_g_f3 l1_Imacec_g pc1,vce(robust) quantile(.`quan')
	eststo Q1M6: quietly qreg Imacec_g_f6 l1_Imacec_g pc1,vce(robust) quantile(.`quan')
	
	eststo Q2M1: quietly qreg Imacec_g pc1, vce(robust) quantile(.`quan')
	eststo Q2M2: quietly qreg Imacec_g l1_Imacec_g pc1,vce(robust) quantile(.`quan')
	eststo Q2M3: quietly qreg Imacec_g_f1 Imacec_g pc1,vce(robust) quantile(.`quan')
	eststo Q2M4: quietly qreg Imacec_g_f2 Imacec_g_f1 pc1,vce(robust) quantile(.`quan')
	eststo Q2M5: quietly qreg Imacec_g_f3 Imacec_g_f2 pc1,vce(robust) quantile(.`quan')
	eststo Q2M6: quietly qreg Imacec_g_f6 Imacec_g_f5 pc1,vce(robust) quantile(.`quan')
	
	esttab Q1M1 Q1M2 Q1M3 Q1M4 Q1M5 Q1M6 using "${tabs}QReg_`quan'_MOMDes_Haztius.tex", replace nogaps b(%-9.3f) se(%-9.3f)  starlevel("*" 0.10 "**" 0.05 "***" 0.01) stats(N, fmt(%12.0gc)  label("Obs.")) label title("Quantile Regressions at `quan'%") nonotes nonumbers style(tex) mtitles("ImacecG t" "ImacecG t" "ImacecG t+1" "ImacecG t+2" "ImacecG t+3" "ImacecG t+6") noline

	esttab Q2M1 Q2M2 Q2M3 Q2M4 Q2M5 Q2M6 using "${tabs}QReg_`quan'__MOMDest-1_Haztius.tex", replace nogaps b(%-9.3f) se(%-9.3f)  starlevel("*" 0.10 "**" 0.05 "***" 0.01) stats(N, fmt(%12.0gc)  label("Obs.")) label title("Quantile Regressions at `quan'%") nonotes nonumbers style(tex) mtitles("ImacecG t" "ImacecG t" "ImacecG t+1" "ImacecG t+2" "ImacecG t+3" "ImacecG t+6") noline

}
drop _est*

* Some preliminary plots
quantile Imacec_g
graph export "${figs}Imacec_g_quartiles.eps",replace

* Quantile regressions for the whole distribution
quietly qreg Imacec_g pc1
grqreg, cons ci ols olsci note("Note: 1PC. This figure plots the Quantile Reg. coef. for each quantile with its respective CI." "The dashed line are OLS coef.")
* Apparently, coeff aren't so different between quantiles.
graph export "${figs}QReg_Coeff_Haztius.eps",replace

* Quantile regressions for the whole distribution, including lag
quietly qreg Imacec_g l1_Imacec_g pc1
grqreg, cons ci ols olsci note("Note: 1PC. This figure plots the Quantile Reg. coef. for each quantile with its respective CI." "The dashed line are OLS coef.")
* Apparently, coeff aren't so different between quantiles.
graph export "${figs}QReg_Coeff_LagImacecG_MOMDes_Haztius.eps",replace

* --------------------------------------------------------------------------- *
* Out of sample predictions
* --------------------------------------------------------------------------- *
gen aux_date = [_n]
* Do a loop for each specification of t+h. Local outside the loop helps to 
* check the counterfactual (AR(1)). The local will set the dep and indep 
* variables since we have two diff specif, with and without lagged econ 
* activity in regressors

* COUNTERFACTUAL
local indep_var l1_Imacec_g
local matrix_name mean_errors_CF
* ESTIMATION WITH FCI
*local indep_var l1_Imacec_g pc1
*local matrix_name mean_errors
local t_h = 0
count
* we begin the loop after half of the sample
local beg_loop = (`r(N)'-4)/2
local end_loop = (`r(N)'-3)
matrix `matrix_name'= J(5,2,.)
foreach depvar of varlist Imacec_g F.Imacec_g F2.Imacec_g F3.Imacec_g F6.Imacec_g {
	gen ofs_pred_ = .
 		forvalues i= `beg_loop'/`end_loop'{
			quietly reg `depvar' `indep_var' if aux_date <= `i'
*			quietly reg `depvar' L.`depvar' `indep_var' if aux_date <= `i'
			local date_pr =`i'+2
			quietly predict prdn if aux_date==`date_pr'
			quietly replace ofs_pred_ = prdn if aux_date==`date_pr'
			drop prdn
		}		
	rename ofs_pred_ ofs_pred_`t_h'
	gen abs_dif_`t_h' = abs(`dep_var' - ofs_pred_`t_h')
	quietly sum abs_dif_`t_h'
	matrix  `matrix_name'[(`t_h'+1),1] = `r(mean)'
	quietly sum abs_dif_`t_h' if aux_date<=357
	matrix `matrix_name'[(`t_h'+1),2] = `r(mean)'
	local t_h = `t_h'+1				
}

mat list `matrix_name'
drop ofs_pred* abs_dif*
* Export results
frmttable using "${tabs}MeanOFSErrors_MOM_Haztius_t2.tex", statmat(mean_errors) sfmt(g) sdec(3) ctitle("Model" "Whole Sample" "Pre-2019m10") rtitle(ImacecG t \ ImacecG t+1\ ImacecG t+2 \ ImacecG t+3 \ ImacecG t+6) fragment tex  hline(11000001) vlines(01000) replace

frmttable using "${tabs}MeanOFSErrors_MOM_Haztius_t2.tex", statmat(mean_errors_CF) sfmt(g) sdec(3) ctitle("Model" "Whole Sample CF" "Pre-2019m10 CF") rtitle(ImacecG t \ ImacecG t+1\ ImacecG t+2 \ ImacecG t+3 \ ImacecG t+6) fragment tex  hline(01000001) vlines(000100) merge

* -------------------------------

* QUANTILE REGRESSIONS OUT OF SAMPLE PREDICTION
* Same idea as before, but with quantile regressions (less stat power)
* I do this because qreg doesn't allow time-series operators, it makes it easy for me to do the loop
*rename l1_Imacec_g Imacec_g_f0
*gen Imacec_g_f0 = Imacec_g
*rename Imacec_g_f5 Imacec_g_f4

* COUNTERFACTUAL
local matrix_name mean_errors_CF
local indep_var l1_Imacec_g
* ESTIMATION WITH FCI
local matrix_name mean_errors
local indep_var l1_Imacec_g pc1
local t_h = 0
count
local beg_loop = (`r(N)'-4)/2
local end_loop = (`r(N)'-3)
matrix `matrix_name'= J(5,2,.)
foreach depvar of varlist Imacec_g Imacec_g_f1 Imacec_g_f2 Imacec_g_f3 Imacec_g_f6 {
	local t_h_1 = `t_h'-1
	quietly gen ofs_pred_ = .
 		forvalues i= `beg_loop'/`end_loop'{
			quietly qreg `depvar' `indep_var' if aux_date <= `i',  vce(robust) quantile(.2)
			local date_pr =`i'+3
			quietly predict prdn if aux_date==`date_pr'
			quietly replace ofs_pred_ = prdn if aux_date==`date_pr'
			drop prdn
		}		
	rename ofs_pred_ ofs_pred_`t_h'
	quietly gen dif_`t_h' = `depvar' - ofs_pred_`t_h'
	quietly gen dummy_dif_`t_h' = 1 if dif_`t_h'<0
	quietly replace dummy_dif_`t_h' = . if dif_`t_h'==.
	quietly replace dummy_dif_`t_h' = 0 if dif_`t_h'>0
	quietly sum dummy_dif_`t_h'
	matrix  `matrix_name'[(`t_h'+1),1] = `r(mean)'
	quietly sum dummy_dif_`t_h' if aux_date<=357
	matrix `matrix_name'[(`t_h'+1),2] = `r(mean)'
	local t_h = `t_h'+1				
}

mat list `matrix_name'
drop ofs_pred* dif* dummy_dif*

frmttable using "${tabs}OFSErrors_QReg_MOMDes_Hatzius.tex", statmat(mean_errors) sfmt(g) sdec(3) ctitle("Model" "Whole Sample" "Pre-2019m10") rtitle(ImacecG t \ ImacecG t+1\ ImacecG t+2 \ ImacecG t+3 \ ImacecG t+6) fragment tex  hline(11000001) vlines(01000) replace

frmttable using "${tabs}OFSErrors_QReg_MOMDes_Hatzius.tex", statmat(mean_errors_CF) sfmt(g) sdec(3) ctitle("Model" "Whole Sample CF" "Pre-2019m10 CF") rtitle(ImacecG t \ ImacecG t+1\ ImacecG t+2 \ ImacecG t+3 \ ImacecG t+6) fragment tex  hline(01000001) vlines(000100) merge
