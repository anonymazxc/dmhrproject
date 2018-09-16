****FIRST STEP: Describe the variables****

count
*There are 4215 rows

sum randid
*There is one row per subject

hist glucose
sum glucose, detail 
*The plot seems slightly skewed to the right but the similar mean and median suggest it follows a normal distribution  
*Exposure has missing data (3,828 vs 4215)

**Categorise glucose by clinically meaningful levels
sum glucose, detail
gen glucat = glucose
recode glucat (39/71.999=1) (72/108=2) (108.1/180=3) (180.1/395=4)
tab glucat
sum glucat, detail
*Low, normal, borderline and high glucose levels (mg/dL)

tab mi_fchd
*602 subjects were hospitalised for MI or fatal CHD during follow-up period

tab death
*1388 deaths during follow-up period (including those for CHD)

tab sex
*There are slightly more men than women

hist bmi 
sum bmi, detail  
*Some missing data  4,198 vs 4215
*The plot seems slightly skewed to the right but the similar mean and median suggest it follows a normal distribution

**Categorise BMI by clinically meaningful levels
sum bmi, detail
gen bmicat = bmi
recode bmicat (15/18.499=1) (18.5/24.999=2) (25/29.999=3) (30/58=4)
tab bmicat
*Underweight, Normal, Overweight, Obese

sum bmicat, detail

**Create variable age
gen age = (dob/-365.25)
hist age
sum age, detail
hist age
*Age ranges from 36 to 74 years

**Categorise Age by quantiles
xtile agecat = age, nquantiles(4)
tab agecat, nol
sum agecat, detail

**Categorise cigarettes per day (Non-smoker=0, Moderate smoker<20, Heavy smoker <70)
sum cigpday, detail
hist cigpday

*Cigpday does not follow a normal distribution. It needs to be transformed

gen cigcat = cigpday
recode cigcat (0=1) (1/20=2) (21/70=3)  
hist cigcat
tab cigcat, nol

*The categories are: "Non-smoker", "Moderate smoker" (1 a 20 per day), "Heavy-smoker" (21 a 70 per day)

sum cigcat, detail

**Describe education and prior hypertension
tab educ
tab prior_hyp

*Describe follow-up time
stset endfollowup, enter(baseline_visit) origin(dob) fail(mi_fchd) id(randid) scale(365.25)
scatter death last_followup
*Follow-up times are different 

tab last_followup if death==0

*Follow up amongst those who did not die during follow up period ended at some point within two years




****DEAL WITH MISSING DATA****

* Count glucose missing values 
codebook glucose

* Create missingness indicator
gen miss_glucose = missing(glucose)
tab miss_glucose 

* Reasons for missingness
logistic miss_glucose i.sex i.agecat i.cigcat i.educ i.prior_hyp i.bmicat, base

* Count BMI missing values 
codebook bmi

* Create missingness indicator
gen miss_bmi = missing(bmi)
tab miss_bmi 

* Reasons for missingness
logistic miss_bmi i.sex i.agecat i.cigcat i.educ i.prior_hyp i.glucat, base

*Variables related to the variables (without missing data)

regress glucose i.sex i.agecat i.cigcat i.educ i.prior_hyp i.bmicat, base

regress bmi i.sex i.agecat i.cigcat i.educ i.prior_hyp i.glucat, base

*Complete records

stcox glucose i.sex i.cigcat, base
stcox, nohr

stcox bmi i.glucat, base
stcox, nohr

*Missing indicator

gen glucose_ind = glucose
replace glucose_ind = 0 if glucose==.
stcox glucose_ind miss_glucose i.sex i.cigcat, base
stcox, nohr

gen bmi_ind = bmi
replace bmi_ind = 0 if bmi==.
stcox bmi_ind miss_bmi i.glucat, base
stcox, nohr

*Implement mean imputation 

gen glucose_imp = glucose
summ glucose
replace glucose_imp = r(mean) if glucose==.
stcox glucose_imp i.sex i.cigcat, base
stcox, nohr

gen bmi_imp = bmi
summ bmi
replace bmi_imp = r(mean) if bmi==.
stcox bmi_imp i.glucat, base
stcox, nohr

*Implement regression imputation

gen glucose_regressimp = glucose
regress glucose i.agecat i.cigcat i.bmicat i.prior_hyp i.sex i.educ
predict glucosepredict
replace glucose_regressimp = glucosepredict if glucose==.
stcox glucose_regressimp i.sex i.cigcat, base
stcox, nohr

gen bmi_regressimp = bmi
regress bmicat i.agecat i.cigcat i.glucat i.prior_hyp i.sex i.educ
predict bmipredict
replace bmi_regressimp = bmipredict if bmi==.
stcox bmi_regressimp i.sex i.cigcat, base
stcox, nohr

*Implement regression imputation including the outcome

gen glucose_regressimpout = glucose
regress glucose i.agecat i.cigcat i.bmicat i.prior_hyp i.sex i.educ i.mi_fchd
predict glucosepredict2
replace glucose_regressimpout = glucosepredict2 if glucose==.
stcox glucose_regressimpout i.sex i.cigcat, base
stcox, nohr

gen bmi_regressimpout = bmi
regress bmicat i.agecat i.cigcat i.glucat i.prior_hyp i.sex i.educ i.mi_fchd
predict bmipredict2
replace bmi_regressimpout = bmipredict2 if bmi==.
stcox bmi_regressimpout i.sex i.cigcat, base
stcox, nohr



****SECOND STEP: CHOOSE AND CREATE A MODEL****


**Is glucose associated with hospitalisation for MI or death for CHD? 

egen endfollowup = rowmin(date_mifchd last_followup)
format endfollowup %td
label var endfollowup "End of follow-up"

strate glucat, per(1000) graph
stmh glucat

*Yes (p<0.001)

*Kaplan-Meier estimates

sts graph, by(glucat)

*Unadjusted cox model

stcox i.glucat

**Are there any confounders or effect modifiers?

*What variables are associated with the outcome?

strate agecat, per(1000) graph
stmh agecat

*Yes (p<0.001)

strate sex, per(1000) graph
stmh sex
*Yes (p<0.001)

strate bmicat, per(1000) graph
stmh bmicat
*Yes (p<0.001)

strate cigcat, per(1000) graph
stmh cigcat
*Yes (p<0.001)

strate educ, per(1000) graph
stmh educ
*No (p=0.2512)

strate prior_hyp, per(1000) graph
stmh prior_hyp
*Yes (p<0.001)

*What variables are also associated with the exposure?

tab glucat agecat, col chi

tab glucat sex, col chi

tab glucat bmicat, col chi

tab glucat cigcat, col chi

tab glucat prior_hyp, col chi

*All of them (p<0.001) except sex (p=0.06).

*Is there evidence for effect modification?

stmh glucose, by(agecat) 
*p=0.17
stmh glucat, by(sex)
*p=0.72
stmh glucat, by(bmicat)
*p=0.55
stmh glucat, by(cigcat)
*p=0.14
stmh glucat, by(prior_hyp)
*p=0.13

*No evidence for effect modification


*********
*Create a model including the confounders and effect modifiers

***Cox model 

stcox i.glucat i.agecat i.sex i.bmicat i.cigcat i.prior_hyp, base
stcox, nohr

*Better model
stcox i.glucat i.sex i.cigcat i.prior_hyp, base
stcox,nohr

*Proportional hazards assumption

estat phtest, detail 

*Non-proportional hazards assumption

stcox i.glucat i.cigcat i.prior_hyp, base strata (sex cigcat)
stcox,nohr

estat phtest, detail 

*Assumption is met

**********


***Explore the nature of the model

**Compare categorical and linear models

findit partpred

stcox glucose i.cigcat i.prior_hyp, base strata (sex)

partpred hr, for(glucose) ref(glucose 78) ci(hr_cl hr_cu) eform

twoway (rarea hr_cl hr_cu glucose, sort pstyle(ci)) ///
	(line hr glucose, sort) ///
	, legend(off) xtitle("Continuous glucose") ytitle("Hazard Ratio")
	
	
stcox ib2.glucat i.cigcat i.prior_hyp, base strata (sex)

partpred hrcat, for(ib2.glucat) ci(hr_cl_cat hr_cu_cat) eform

twoway (rcap hr_cl_cat hr_cu_cat glucat, sort pstyle(ci)) ///
	(scatter hrcat glucat, sort) ///
	, legend(off) xtitle("Categorised glucose") ytitle("Hazard Ratio")
	
twoway (line hrcat glucose, sort) ///
(line hr glucose, sort) ///
, legend(off) xtitle("Glucose") ytitle("Hazard Ratio")

*The models diverge from values above 230 

stcox glucose i.cigcat i.prior_hyp, base strata (sex)
stcox, nohr
est store a
stcox i.glucat i.cigcat i.prior_hyp, base strata (sex)
stcox, nohr
est store b
lrtest a b 

*The two models are not significantly different (p=1.00)


**Polynomial model (quadratic vs. linear)
 
stcox glucose i.cigcat i.prior_hyp, base strata (sex)
stcox, nohr
est store a
stcox glucose##glucose i.cigcat i.prior_hyp, base strata (sex)
predict mi_fchd_quad
est store b
lrtest a b  

*p<0.001. There is evidence for departures from linearity
*Quadratic is a better model

twoway (scatter mi_fchd glucose)(line mi_fchd_quad glucose)


**Polynomial model (cubic vs. linear)
 
stcox glucose i.cigcat i.prior_hyp, base strata (sex)
est store a
stcox glucose##glucose##glucose i.cigcat i.prior_hyp, base strata (sex)
est store b
lrtest a b     

*p<0.001. There is evidence for departures from linearity

stcox glucose##glucose##glucose i.cigcat i.prior_hyp, base strata (sex)
predict mi_fchd_quad2
twoway (scatter mi_fchd glucose)(line mi_fchd_quad2 glucose)


*Cubic is a better model

*Plot quadratic and cubic
twoway (scatter mi_fchd_quad glucose)(line mi_fchd_quad2 glucose)

*They seem to provide similar estimates
 
**Log transformation

gen log_glucose= ln(glucose)

stcox glucose i.cigcat i.prior_hyp, base strata (sex)
est store a
stcox log_glucose glucose i.cigcat i.prior_hyp, base strata (sex)
est store b
lrtest a b 

*Log is better than linear (p<0.05)

*Is log better than quadratic / cubic?

*Quadratic vs. log

stcox glucose##glucose i.cigcat i.prior_hyp, base strata (sex)
est store a
stcox log_glucose glucose i.cigcat i.prior_hyp, base strata (sex)
est store b
lrtest a b 

*Quadratic is a better model (p<0.001)

*Cubic vs. log

stcox glucose##glucose##glucose i.cigcat i.prior_hyp, base strata (sex)
est store a
stcox log_glucose glucose i.cigcat i.prior_hyp, base strata (sex)
est store b
lrtest a b 

*Cubic is a better model (p<0.001)






