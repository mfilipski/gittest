*date: Aug 24,2014
*goal: psychological effect of the earthquake on consumption behavior
clear
set more off

* local files
local curdir "D:\Docs\ChinaEarthquake\SavingsRate\mywork\dofiles"
local rawdat "D:\Docs\ChinaEarthquake\SavingsRate\mywork\Data\JinlingData"
local madedat "D:\Docs\ChinaEarthquake\SavingsRate\mywork\Data\madedat"
local output "D:\Docs\ChinaEarthquake\SavingsRate\mywork\output"

/*
* Table 1: without the controls:  
*---------------------------------
* tables in my first draft: 
use "`rawdat'\Saving Rate_FD.dta",clear
tab cm sm
count
reg lnsavingrate_0709 lndistance3_2   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_dis3_2_t1 
reg lnsavingrate_0709 lndistance5  ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_dis5_t1 
reg lnsavingrate_0711 lndistance3_2  ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm) 
est store lnsrate0711_dis3_2_t1
reg lnsavingrate_0711 lndistance5  ///
	age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm) 
est store lnsrate0711_dis5_t1

* use estout for temporary output, outreg2 for the excel ones
estout lnsrate0709_dis3_2_t1 lnsrate0709_dis5_t1 lnsrate0711_dis3_2_t1 lnsrate0711_dis5_t1   /// 
	, cell(b(star fmt(%9.3f))) stats(N r2_a) var(9) model(12)  ///
	order(lndistance3_2 lndistance5) drop(age* sex* hh* educ* labor* party* acreag*) 
	
* lntnincomef0709 lntnincomef07 cpi0709 giftexchange0709 loan0709
* lntnincomef0711 lntnincomef07 cpi0711 giftexchange0711 loan0711

* Adding the specifications with damages
*tab cm sm 
*use "`rawdat'\Human Injury and House Damage_Village Level.dta" 
*list 
*crash

merge m:1 sm cm using "`rawdat'\Human Injury and House Damage_Village Level.dta" 
tab _m
sum
tab cm percent_house_collapse

reg lnsavingrate_0709 percent_house_collapse   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
reg lnsavingrate_0709 percent_house_impaired   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
reg lnsavingrate_0709 percent_human_died   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
reg lnsavingrate_0709 percent_human_injured   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)


* Note: This recode is important - makes a difference in results significance, 
* Not sure the assumption is true for human injuries 
recode percent_house*  percent_human* (.=0)

reg lnsavingrate_0709 percent_house_collapse   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_pccol_t1
reg lnsavingrate_0709 percent_house_impaired   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_pcimp_t1
reg lnsavingrate_0709 percent_human_died   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_pcdie_t1
reg lnsavingrate_0709 percent_human_injured   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_pcinj_t1

reg lnsavingrate_0711 percent_house_collapse   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0711_pccol_t1
reg lnsavingrate_0711 percent_house_impaired   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0711_pcimp_t1
reg lnsavingrate_0711 percent_human_died   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0711_pcdie_t1
reg lnsavingrate_0711 percent_human_injured   ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0711_pcinj_t1
estout lnsrate0709_pccol_t1 lnsrate0709_pcimp_t1 lnsrate0709_pcdie_t1 lnsrate0709_pcinj_t1 ///
	, cell(b(star fmt(%9.3f))) stats(N r2_a) var(9) model(12)  ///
	order(percent*) drop(age* sex* hh* educ* labor* party* acreag*) 

estout lnsrate0711_pccol_t1 lnsrate0711_pcimp_t1 lnsrate0711_pcdie_t1 lnsrate0711_pcinj_t1 ///
	, cell(b(star fmt(%9.3f))) stats(N r2_a) var(9) model(12) ///
	order(percent*) drop(age* sex* hh* educ* labor* party* acreag*) 


outreg2 [lnsrate0709_dis3_2_t1 lnsrate0709_dis5_t1 lnsrate0709_pcinj_t1 lnsrate0709_pccol_t1 ///
		lnsrate0711_dis3_2_t1 lnsrate0711_dis5_t1 lnsrate0711_pcinj_t1 lnsrate0711_pccol_t1 ]  /// 
	using "`output'\table1_lnYC_FD", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) ///
	sortvar(lndistance3_2 lndistance5) drop(age* sex hh* educ labor party acreag) nocons


* Table 2: FE specification and Adding controls : 
*---------------------------------------------------------------
use "`rawdat'\Saving Rate_FE.dta",clear
set more off
tset id year,yearly
* Fixed effects no controls: 
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis3_2_t2

* Fixed effects with extra controls:
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	lntnincomef  /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis3_2_t2_c1
xtreg lnincome_consumption   lndistance3_2 distime093_2 distime113_2 /// 
	cpi  /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis3_2_t2_c2
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	giftexchangef loanf /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis3_2_t2_c3
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	lntnincomef cpi giftexchangef loanf /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis3_2_t2_c4

* outreg2 [lnsrate_dis3_2_t2 lnsrate_dis3_2_t2_c1 lnsrate_dis3_2_t2_c2 lnsrate_dis3_2_t2_c3 lnsrate_dis3_2_t2_c4] ///
*	using "`output'\table2_controls_FE", ///
*	adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) ///
*	sortvar(distime093_2 distime095 distime113_2 distime115 lntnincomef cpi giftexchangef loanf )  ///
*	drop(age* sex hh* educ labor party acreag) nocons
	
* Try the same with random effects (=doesn't drop the distance variable)
set more off
tset id year,yearly
* Fixed effects no controls: 
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,re vce(cluster cm)
est store lnsrate_dis3_2r_t2

* Fixed effects with extra controls:
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	lntnincomef  /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,re vce(cluster cm)
est store lnsrate_dis3_2r_t2_c1
xtreg lnincome_consumption   lndistance3_2 distime093_2 distime113_2 /// 
	cpi  /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,re vce(cluster cm)
est store lnsrate_dis3_2r_t2_c2
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	giftexchangef loanf /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,re vce(cluster cm)
est store lnsrate_dis3_2r_t2_c3
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 /// 
	lntnincomef cpi giftexchangef loanf /// 
	age age2 sex educ hhsize hhsize2 labor acreage_land partymember,re vce(cluster cm)
est store lnsrate_dis3_2r_t2_c4

* we take r-squared overall in RE specifications
outreg2 [lnsrate_dis3_2r_t2 lnsrate_dis3_2r_t2_c1 lnsrate_dis3_2r_t2_c2 lnsrate_dis3_2r_t2_c3 lnsrate_dis3_2r_t2_c4] ///
	using "`output'\table2_controls_RE", ///
	adds(Adjusted R-square, e(r2_o)) replace  excel dec(3) ///
	sortvar(distime093_2 distime095 distime113_2 distime115 lntnincomef cpi giftexchangef loanf )  ///
	drop(age* sex hh* educ labor party acreag) nocons

* All together
outreg2 [lnsrate_dis3_2_t2 lnsrate_dis3_2_t2_c1 lnsrate_dis3_2_t2_c2 ///
	lnsrate_dis3_2_t2_c3 lnsrate_dis3_2_t2_c4 lnsrate_dis3_2r_t2_c4] ///
	using "`output'\table2_controls_FERE", ///
	adds(Adjusted R-square, e(r2_o)) replace  excel dec(3) ///
	sortvar(distime093_2 distime095 distime113_2 distime115 lntnincomef cpi giftexchangef loanf lndistance3_2 )  ///
	drop(age* sex hh* educ labor party acreag o.lndistance3_2) nocons

	

* Controlling that incomes and prices are not correlated to distance 
* Method 1: correlation matrix
corr distance3_2 lntnincomef cpi

*/

* Method 2: bar charts over affected-unaffected: 
use "`rawdat'\Saving Rate_FE.dta",clear
merge m:1 cm using `madedat'\vil_coord_magn
tab _m
gen affected = (magn==2 |magn==3)
gen tnincomef = exp(lntnincomef)
graph bar (mean) lntnincome ,over(magn)
graph bar (mean) tnincome ,over(magn)
graph bar (mean) tnincome ,over(magn) over(year)
graph bar (mean) lntnincome ,over(affected) over(year)
cibar lntnincome, over1(affected) over2(year)
cibar lntnincome, over1(magn) over2(year)
graph bar (mean) tnincome ,over(affected)
graph bar (mean) cpi ,over(magn)
graph bar (mean) cpi ,over(affected)

graph twoway (bar (mean) cpi) | (bar (mean) lntn) , over(affected)
*twoway scatter cpi distance3
crash


	
crash 
version 12
//1. first-difference model
//1.1 saving rate
use "`rawdat'\Saving Rate_FD.dta",clear
set more off

reg lnsavingrate_0709 lndistance3_2 lntnincomef0709 lntnincomef07 cpi0709 giftexchange0709 loan0709 ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_dis3_2_ex 
reg lnsavingrate_0709 lndistance5 lntnincomef0709 lntnincomef07 cpi0709 giftexchange0709 loan0709 ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm)
est store lnsrate0709_dis5_ex 

reg lnsavingrate_0711 lndistance3_2 lntnincomef0711 lntnincomef07 cpi0711 giftexchange0711 loan0711 ///
	age07 age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm) 
est store lnsrate0711_dis3_2_ex
reg lnsavingrate_0711 lndistance5 lntnincomef0711 lntnincomef07 cpi0711 giftexchange0711 loan0711 age07 ///
	age07_2 sex07 hhsize07 hhsize07_2 educ07 labor07 partymember07 acreage_land07,vce(cluster cm) 
est store lnsrate0711_dis5_ex

outreg2 [lnsrate0709_dis3_2_ex lnsrate0709_dis5_ex lnsrate0711_dis3_2_ex lnsrate0711_dis5_ex]  /// 
	using "`output'\saving rate_log2_ex", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 


//1.2 playing majiang
use "`rawdat'\Playing Majiang_FD.dta",clear
set more off
reg hh_frequency0711 lndistance3_2 lntnincomef lntnincomef_2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 /// 
	labor07 acreage_land07 partymember07,vce(cluster cm)
est store frequency_dis3_2_ex
reg hh_frequency0711 lndistance5 lntnincomef lntnincomef_2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 /// 
	labor07 acreage_land07 partymember07,vce(cluster cm)
est store frequency_dis5_ex

reg hh_losingamount0711 lndistance3_2 lntnincomef lntnincomef_2 age07 age07_2 sex07 educ07 hhsize07 ///
	hhsize07_2 labor07 acreage_land07 partymember07,vce(cluster cm)
est store losingamount_dis3_2_ex 
reg hh_losingamount0711 lndistance5 lntnincomef lntnincomef_2 age07 age07_2 sex07 educ07 hhsize07 ///
	hhsize07_2 labor07 acreage_land07 partymember07,vce(cluster cm)
est store losingamount_dis5_ex 

outreg2 [frequency_dis3_2_ex frequency_dis5_ex losingamount_dis3_2_ex losingamount_dis5_ex] /// 
	using "`output'\playing majiang_for thesis", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 

//1.3 alcohol and cigarette consumption
use "`rawdat'\Alcohol and Cigarette Consumption_FD.dta",clear
set more off
reg alcohol0709 lndistance3_2 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 /// 
	acreage_land07 partymember07,vce(cluster cm)
est store alcohol0709_dis3_2 
reg alcohol0709 lndistance5 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 /// 
	acreage_land07 partymember07,vce(cluster cm)
est store alcohol0709_dis5 

reg alcohol0711 lndistance3_2 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 /// 
	acreage_land07 partymember07,vce(cluster cm)
est store alcohol0711_dis3_2 
reg alcohol0711 lndistance5 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 /// 
	acreage_land07 partymember07,vce(cluster cm)
est store alcohol0711_dis5 

*alcohol0709_humi alcohol0709_houc alcohol0711_humi alcohol0711_houc
outreg2 [alcohol0709_dis3_2 alcohol0709_dis5 alcohol0711_dis3_2 alcohol0711_dis5 ] /// 
	using "`output'\alcohol", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 

reg cigarette0709 lndistance3_2 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 ///
	acreage_land07 partymember07,vce(cluster cm)
est store cigarette0709_dis3_2 
reg cigarette0709 lndistance5 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 ///
	acreage_land07 partymember07,vce(cluster cm)
est store cigarette0709_dis5 

reg cigarette0711 lndistance3_2 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 ///
	acreage_land07 partymember07,vce(cluster cm)
est store cigarette0711_dis3_2 
reg cigarette0711 lndistance5 lntnincomef lntnincomef2 age07 age07_2 sex07 educ07 hhsize07 hhsize07_2 labor07 ///
	acreage_land07 partymember07,vce(cluster cm)
est store cigarette0711_dis5 

*cigarette0709_dis1 cigarette0709_humi cigarette0709_houc cigarette0711_dis1 cigarette0711_humi cigarette0711_hou
outreg2 [ cigarette0709_dis3_2 cigarette0709_dis5 cigarette0711_dis3_2 cigarette0711_dis5] ///
	using "`output'\cigarette", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 



//1.4 consumption on consumer durables
// only performing fixed-effect model

//2. fixed-effect model
//2.1 saving rate
use "`rawdat'\Saving Rate_FE.dta",clear
set more off
tset id year,yearly
xtreg lnincome_consumption time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef cpi giftexchangef loanf age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis3_2
xtreg lnincome_consumption time09 time11 lndistance5 distime095 distime115 lntnincomef cpi giftexchangef loanf age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store lnsrate_dis5

outreg2 [lnsrate_dis3_2 lnsrate_dis5] using "`output'\saving rate_fe_for submission", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 

//2.2 playing majiang
use "`rawdat'\Playing Majiang_FE.dta",clear
set more off
tset id year,yearly
xtreg hh_frequency time11 lndistance3_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store frequency_dis3_2
xtreg hh_frequency time11 lndistance5 distime115 lntnincomef lntnincomef2  age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store frequency_dis5

xtreg hh_losingamountf time11 lndistance3_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store losingamount_dis3_2
xtreg hh_losingamountf time11 lndistance5 distime115 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store losingamount_dis5

outreg2 [frequency_dis3_2 frequency_dis5 losingamount_dis3_2 losingamount_dis5] using "`output'\playing majiang_fe", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 

//2.3 alcohol and cigarette consumption
use "`rawdat'\Alcohol and Cigarette Consumption_FE.dta",clear
set more off
tset id year,yearly

xtreg lnalcohol time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store alcohol_dis3_2
xtreg lnalcohol time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store alcohol_dis5

xtreg lncigarette time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store cigarette_dis3_2
xtreg lncigarette time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store cigarette_dis5

outreg2 [alcohol_dis3_2 alcohol_dis5 cigarette_dis3_2 cigarette_dis5] using "`output'\alcohol cigarette_fe", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 

//2.4 consumption on consumer durables
//a. yearly expenditure on consumer durables
use "`rawdat'\Yearly Expenditure on Consumer Durables_FE.dta", clear
set more off
tset id year,yearly

xtreg lncdyrexpf time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store cdexp_32
xtreg lncdyrexpf time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store cdexp_5

outreg2 [cdexp_32 cdexp_5] using "`output'\consumer durable_expenditure", adds(Adjusted R-square, e(r2_a)) replace  excel dec(3) 


//units of consumer durables owned
use "`rawdat'\Units of Consumer Durables Owned_FE.dta",clear
set more off
tset id year,yearly

gen lncondurable=ln(condurable+1)
xtreg lncondurable time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store condurable_dis3_2
xtreg lncondurable time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(cluster cm)
est store condurable_dis5  
outreg2 [ condurable_dis3_2 condurable_dis5] using "`output'\consumer durable_aggregate", adds( Adjusted R-square, e(r2_a) ) replace  excel dec(3) 

xtpoisson condurable time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store condurable_32_poisson
xtpoisson condurable time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store condurable_5_poisson                
* don't replace so that it appends:              
outreg2 [ condurable_32_poisson condurable_5_poisson] using "`output'\consumer durable_aggregate", adds( chi-square test, e(chi2) )  excel dec(3) 
*condurable_dis3_2 condurable_dis5 ,

*esttab  condurable_32_poisson condurable_5_poisson using "`output'\consumer durable_aggregatetab.csv", stats( chi2 ) replace  


//consumer durables for entertainment
gen cd_enter=tv+cassette+stereo+camero+vcr+videocamero+dvdplayer+computer
gen cd_trans=bicycle+motocycle+car+phone+cell
gen cd_comf=washingma+sewingma+refri+microwave+elecooker+furniture+aircon+waterheater+fan

xtpoisson cd_enter time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex ///
	educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store cd_enter32_poisson
xtpoisson cd_enter time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex ///
	educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store cd_enter5_poisson                

xtpoisson cd_trans time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex ///
	educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store cd_trans32_poisson
xtpoisson cd_trans time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex /// 
	educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store cd_trans5_poisson                

xtpoisson cd_comf time09 time11 lndistance3_2 distime093_2 distime113_2 lntnincomef lntnincomef2 age age2 sex ///
	educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store cd_comf32_poisson
xtpoisson cd_comf time09 time11 lndistance5 distime095 distime115 lntnincomef lntnincomef2 age age2 sex /// 
	educ hhsize hhsize2 labor acreage_land partymember,fe vce(robust)
est store cd_comf5_poisson                

outreg2 [cd_enter32_poisson cd_enter5_poisson cd_trans32_poisson cd_trans5_poisson cd_comf32_poisson cd_comf5_poisson] ///
	using "`output'\consumer durable_aggregate 2", adds(chi-square test, e(chi2) ) replace  excel dec(3) 






