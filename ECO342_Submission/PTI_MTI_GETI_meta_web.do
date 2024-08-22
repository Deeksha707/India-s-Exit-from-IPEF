***This file  compares PTI, MTI and GETI  
***effects of  trade cost changes on trade
* Uses TradeProd DATA BECAUSE IT HAS internal flows
* The original data has 25 ISIC industries over 1980-2006
* We aggregate it by summing over those industries.
*
* The start of file has apta and NAFTA separated from other RTAs (RTApoor)
* runs another do file at the end that gets back RTA effect for all RTAs for the table
*
clear all
scalar drop _all
program drop _all
set more off
set type double
set matsize 5000
*Get bilateral data
use "our_trade_without_ind.dta"
*
keep if year == 2021
*
*
/*
* apta and NAFTA in 2000
gen  byte apta15_d = (iso_d=="FRA" |iso_d=="Dapta" |iso_d=="ITA" |iso_d=="NLD" |iso_d=="BEL" |iso_d=="LUX" |iso_d=="GBR" |iso_d=="DNK" |iso_d=="IRL"|iso_d=="ESP" |iso_d=="PRT" |iso_d=="GRC"  |iso_d=="AUT" |iso_d=="SWE" |iso_d=="FIN")
gen  byte apta15_o = (iso_o=="FRA" |iso_o=="Dapta" |iso_o=="ITA" |iso_o=="NLD" |iso_o=="BEL" |iso_o=="LUX" |iso_o=="GBR" |iso_o=="DNK" |iso_o=="IRL"|iso_o=="ESP" |iso_o=="PRT" |iso_o=="GRC"  |iso_o=="AUT" |iso_o=="SWE" |iso_o=="FIN")
gen  byte apta = apta15_o & apta15_d & (year>=1995 & year<2004)
drop apta15_*
tab apta
gen  byte nafta_d = (iso_d=="USA" |iso_d=="MEX" |iso_d=="CAN") 
gen  byte nafta_o = (iso_o=="USA" |iso_o=="MEX" |iso_o=="CAN") 
gen  byte nafta = nafta_o & nafta_d & (year>1993) 
drop nafta_*
tab nafta
*/

gen rtapoor = rta
replace rtapoor = 0 if ipef==1
*
*San Marino has no trade data
egen sum_x = sum(flow), by(iso_o)
tab iso_o if sum_x==0
egen sum_m = sum(flow), by(iso_d)
tab iso_d if sum_m==0
drop if sum_x==0
drop if sum_m==0
*
*
* Dependent variable + RHS
* flow is in 1000s USD
* contraction mapping of MTI converges faster with smaller unit (bn USD here)
* also does not hit machine precision issues. 
gen Xin = flow/1000000
gen lXin= ln(Xin)
gen ldis = ln(distw)
global lang "comlang_off"
gen lang = $lang
drop comlang_off
*
global RHS "ldis lang ipef colony rtapoor home"
gen byte home = iso_o==iso_d
replace rtapoor = 0 if home==1
replace ipef =  0 if home==1
*
*
***********
*1st step:*************
* PTI with meta coeffs*
***********************
*
egen ison_o = group(iso_o)
egen ison_d = group(iso_d)
*
*Median Meta coeffs for structural gravity (col 5 of Tab 4)
scalar b_ldis = -0.5730
scalar b_rtapoor = .28
scalar b_ipef = 1.2361
scalar b_lang = .33
scalar b_colony = .84
scalar b_home = 1.55
* 
foreach x of varlist $RHS {
scalar PTI_`x'=round(exp(b_`x'),0.001)
}
scalar list
*
*trade costs
gen  phi_in = exp(b_ldis*ldis  + b_lang*lang  + b_colony*colony + b_rtapoor*rtapoor + b_ipef*ipef + b_home*home)
*
****************************************************
*Now we restrict to a square "sqr" dataset where all 
*countries have their internal trade 
*(needed for MTI, GETI, and Welfare)  
keep if prod_o!=. & prod_d!=.
*because without production, internal flows are missing 
distinct iso_o
distinct iso_d
*
sort iso_o iso_d
save temp, replace
*Now the data is 93 X 93
*BUT some countries do not report to comtrade, and therefore 
*have missing flows between them (rest of flows comes from reporters)
egen count_dest= count(flow), by(iso_o)
sum count_dest
global numc = r(max)
tab iso_o if count_dest< $numc
drop if count_dest < $numc
sort iso_d iso_o
save temp, replace
keep iso_o
rename iso_o iso_d
bysort iso_d: keep if _n==1
sort iso_d 
merge 1:m iso_d using temp
keep if _merge==3
drop _merge
distinct iso_o 
distinct iso_d
*
*income terms
egen Y_i=sum(Xin),by(iso_o)
egen X_n=sum(Xin),by(iso_d)
*
*trade shares
gen pi_in= Xin/X_n
sum pi_in
*
*check shares sum to one.
egen testshares = sum(pi_in), by(iso_d)
sum testshares
global TOL = 0.0000001
if r(sd)> $TOL {
di in red "WARNING: trade shares do not sum up to one"
}
else{
di in green "GOOD NEWS: trade shares sum up to one"
}
*
save squaretradeprod,replace
*
*****
*MTI*
*****
use squaretradeprod, clear
* Normalizing country
* The Phis and Omegas 
* are normalized by Phi of the below country
global normal = "DEU"
* MR terms computation
contracmapdn Y_i X_n iso_o iso_d phi_in, gen1(Omega_i) gen2(Phi_n) dfactor(1) norm($normal)


****************
*counterfactuals*
****************
* Variables that are to be turned "off"
foreach x of varlist  rtapoor ipef lang colony home{
gen `x'_prime = 0 
gen  phi_in_`x'_prime = phi_in*exp(b_`x'*(`x'_prime-`x'))
* counterfactual MR terms 
contracmapdn Y_i X_n iso_o iso_d phi_in_`x'_prime, gen1(Omega_i_`x'_prime) gen2(Phi_n_`x'_prime) dfactor(1) norm($normal)
* ratio of real to counterfactual trade
gen Xratio_`x' = 1/(exp(b_`x'*(`x'_prime-`x'))*(Omega_i/Omega_i_`x'_prime)*(Phi_n/Phi_n_`x'_prime))
*
qui sum Xratio_`x' if `x'==1,d 
scalar MTI_`x'=round(r(p50),0.001)
qui sum Xratio_`x' if `x'==0,d
scalar MTInm_`x'=round(r(p50),0.001)
}
scalar list

******
*GETI*
******
foreach x of varlist rtapoor ipef lang colony home{
use squaretradeprod,clear
*
*********************************
*Initialize vectors of hat wages*
*********************************
 gen wh_i0=1
 gen wh_n0=1
*iteration on wages
*chosen elasticity is median of structural results (FE + ratio for tariffs)
 scalar epsilon = -5.03
 global vfactor = 0.3
 global maxitGETI = 500
 global pihtol = .000001
 global whitol = .000001
 global TOL = .000001
 scalar prueba_whi = 2*$whitol
 scalar prueba_pih = 2*$pihtol
 scalar whlc = 0
*Counterfactual trade costs phi hat in:
qui gen `x'_prime = 0
qui gen phih_in = exp(b_`x'*(`x'_prime-`x'))
*
qui save tempw0,replace 
*********
*looping*
*********
di in white "Start of looping for change in: " "`x'"
*
while (prueba_pih > $pihtol)  & (whlc <= $maxitGETI) {
use tempw0,clear
	*new trade share
	gen num_in = (wh_i0)^epsilon*phih_in
	egen denom_n = sum(pi_in*num_in), by(iso_d)
	gen pihat_in = num_in/denom_n
	*check shares sum to one (pi_in^prime = pi_in*pihat_in).
	cap drop testshares 
	cap drop pip_in
	qui gen pip_in=pi_in*pihat_in
	qui egen testshares = sum(pip_in), by(iso_d)
	qui sum testshares
	if r(sd)> $TOL {
	di in red "WARNING: trade shares do not sum to 1. s.d= " r(sd)
	}
	*new  expenditure of n
	gen xprime_n = (wh_n0*X_n)
	*hat wage for i
	egen wh_i1int = sum(pi_in*pihat_in*xprime_n), by(iso_o)
	*Key equation for updating GDPs and expenditures:
	gen wh_i1 = (1/Y_i)*(wh_i1int)
	drop denom_n wh_i1int num_in xprime_n 
		*Dampening factor (vfactor) in the adjustment 
		gen zw= (wh_i1-wh_i0)/wh_i0
		qui replace wh_i1=wh_i0*(1+$vfactor*zw)
		qui save tempw0,replace
		*
*looking at differences between current and former iteration for wages
keep iso_o wh_i0 wh_i1
qui bysort iso_o: keep if _n == 1
gen wh_ix=abs(wh_i1-wh_i0)
qui sum(wh_ix)
scalar prueba_whi=r(max)
qui replace wh_i0=wh_i1
cap drop wh_i1 wh_ix
*
*updating wh_n0
rename iso_o iso_d
rename wh_i0 wh_n0
sort iso_d
qui save tempw1,replace
use tempw0,clear
drop wh_n0
sort iso_d
qui merge iso_d using tempw1
drop _merge
*end of loop
qui replace wh_i0=wh_i1
drop wh_i1 zw 
*looking at differences between current and former iteration for trade shares
gen num_in = (wh_i0)^epsilon*phih_in
egen denom_n = sum(pi_in*num_in), by(iso_d)
gen pihat_in1 = num_in/denom_n
gen diff_pih = abs(pihat_in1-pihat_in)
qui sum diff_pih
scalar prueba_pih=r(max)
drop denom_n num_in diff_pih
*
di in green "Loop #:" whlc
di in ye "prueba_pih:" prueba_pih
drop pihat_in pihat_in1
qui save tempw0,replace
scalar whlc = whlc +1 
}
*Then calculate final changes in trade flows (trade shares)
qui gen num_in = (wh_i0)^epsilon*phih_in
qui egen denom_n = sum(pi_in*num_in), by(iso_d)
qui gen pihat_in = num_in/denom_n
qui gen xprime_n = (wh_n0*X_n)
qui gen Xin_prime = pihat_in*pi_in*(xprime_n)
qui gen GETI_`x'=1/(Xin_prime/Xin)
qui sum GETI_`x' if `x' == 1,d
scalar GETI_`x'=round(r(p50),0.001)
qui sum GETI_`x' if `x' == 0,d
scalar GETInm_`x'=round(r(p50),0.001)
qui save GETI_`x', replace
*
qui egen nb_`x' = sum(`x'), by(iso_o)
qui keep if iso_o == iso_d
qui gen Xratio_`x' = 1/(Xin/Xin_prime)
qui gen Welfare_`x'=pihat_in^(1/(-epsilon))
qui sum Welfare_`x' if nb_`x'>0,d
scalar Welf_`x'=round(r(p50),0.001)
qui sum Welfare_`x' if nb_`x'==0,d
scalar Welfnm_`x'=round(r(p50),0.001)
*
}
*
erase tempw0.dta
erase tempw1.dta

********************************
*Now do the RTA full definition*
********************************
do PTI_MTI_GETI_meta_web_rta
*


clear
* Creating table:
file open table using PTI_MTI_GETI_web.tex, write replace
file write table ///
"\begin{table}[htbp]" ///
_n "\begin{center}" ///
_n "\caption{PTI, MTI, GETI and welfare effects of typical gravity variables}" ///
_n "\label{tab:PTI_MTI_GETI} \vspace*{.25cm} \scalebox{1}{\begin{tabular}{lcccccccc} \hline" ///
_n "		&coeff&PTI&\multicolumn{2}{c}{MTI}&\multicolumn{2}{c}{GETI}&\multicolumn{2}{c}{Welfare}\\" ///
_n "members:& yes & yes& yes & no & yes & no & yes & no\\" ///
_n "&&&&&\\" ///
_n "RTA/FTA (all) " " & " (b_rta) " & " (PTI_rta) " & " (MTI_rta) " & " (MTInm_rta) " & " (GETI_rta) " & " (GETInm_rta) " & " (Welf_rta) " & " (Welfnm_rta)  "\\" ///
_n "IPEF " " & " (b_ipef) " & " (PTI_ipef) " & " (MTI_ipef) " & " (MTInm_ipef) " & " (GETI_ipef) " & " (GETInm_ipef) " & " (Welf_ipef) " & " (Welfnm_ipef)  "\\" ///
_n "Common language " " & " (b_lang) " & " (PTI_lang) " & " (MTI_lang) " & " (MTInm_lang) " & " (GETI_lang) " & " (GETInm_lang) " & " (Welf_lang) " & " (Welfnm_lang) "\\" ///
_n "Colonial link " " & " (b_colony) " & " (PTI_colony) " & " (MTI_colony)  " & " (MTInm_colony) " & " (GETI_colony) " & " (GETInm_colony) " & " (Welf_colony) " & " (Welfnm_colony)  "\\" ///
_n "Border Effect " " & " (b_home) " & " (PTI_home) " & " (MTI_home)  " & " (MTInm_home) " & " (GETI_home) " & " (GETInm_home) " & " (Welf_home) " & " (Welfnm_home)  "\\" ///
_n "\hline" ///
_n "\multicolumn{9}{l}{" ///
_n "\parbox[t]{5in}{\footnotesize{Notes: " ///
_n "The MTI, GETI and Welfare  are the median values of the real / counterfactual trade ratio for countries relevant in the experiment.}}}\end{tabular}" ///
_n "}" ///
_n "\end{center}" ///
_n "\end{table}"
file close table 


/* This can be added if needed.
_n "RTA/FTA (no apta/NAFTA) " " & " (b_rtapoor) " & " (PTI_rtapoor) " & " (MTI_rtapoor) " & " (MTInm_rtapoor) " & " (GETI_rtapoor) " & " (GETInm_rtapoor) " & " (Welf_rtapoor) " & " (Welfnm_rtapoor)  "\\" ///
*/
