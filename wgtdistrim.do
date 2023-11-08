*********************************************************************************
*******Trimming of weights: weight distribution approach (Potter 1990)***********
*********************************************************************************
*** Autor: 				Sebastian Lang                						  ***
*** Mail: 				lang@dzhw.eu 										  ***
*** Mail: 				contact@sebastianlang.eu							  ***
*********************************************************************************
*********************************************************************************

/*
**# 	Fake-weights for testing
clear all
set obs 100000
gen wgt = runiform()*100

local 		wgtvar = `"wgt_kal2"' //var with weights to be trimmed
local 		cutpoint = .001 //upper and lower portion of distribution that should be trimmed (symmetric)
local 		rounds = 10 //iterations
*/

**#		Collect information in the trimming process
mat 		trimdok = J(`rounds'+1,4,0)
mat 		colnames trimdok = min max mean var

**# 	set upper and lower thresholds (at the moment symmetrically)
local 		uppercut = 1-`cutpoint'
local 		lowercut = `cutpoint'

**#		Trimming procedure iterated `rounds' times
forvalue 	round = 1(1)`rounds' {
	**#		Prepare for trimming
	if 			"`round'" == "1" {
		* 	Normalize weights on samplesize (sum of weights == samplesize)
		cap drop 	`wgtvar'_norm
		cap drop 	`wgtvar'_trim
		qui sum 	`wgtvar'
		gen		 	`wgtvar'_norm = `wgtvar' / r(mean)
		clonevar 	`wgtvar'_trim = `wgtvar'_norm

		*		Collect information on normalizes, untrimmed weights
		qui sum 	`wgtvar'_norm
		mat 		trimdok[`round',1] = r(min)
		mat 		trimdok[`round',2] = r(max)
		mat 		trimdok[`round',3] = r(mean)
		mat 		trimdok[`round',4] = r(Var)
	}
	else 		{
		replace 	`wgtvar'_norm = `wgtvar'_trim
	}
	
	*	Clear mata
	mata: 		mata clear
	
	*pause
	
	*	Prevent including obs in N in case of missings on the wgtvar
	preserve
	drop 		if missing(`wgtvar')
	
	**#		Calculation of the shape parameters for the beta distribution
	tempvar		N
	qui sum 	`wgtvar'
	gen 		`N' = r(N)
	putmata 	wi = `wgtvar'_norm, omit replace
	putmata 	N = `N', omit replace
	mata: 		N = N[1,1]
	*disp 		"N"
	*mata:		N
	mata: 		I = J(rows(wi),1,1)
	*disp 		"I"
	*mata:		I
	mata: 		w_quer = (wi'*I) :/ N
	*disp 		"w_quer"
	*mata: 		w_quer
	mata: 		s_sq = (((wi :- w_quer) :* (wi :- w_quer))' * I) / N
	*disp 		"s_sq"
	*mata: 		s_sq
	*V1 - following Potter (1990)
	mata: 		alpha = ((w_quer * ((N * w_quer) - 1)) / (N * s_sq)) + 2
	mata: 		beta = ((N * w_quer) - 1) * (((w_quer * ((N * w_quer) - 1)) / (N * s_sq)) + 1)
	*V2 - more general formulas for calcuating alpha and beta lead to strange results
	*mata: 		alpha = (w_quer * (s_sq + (w_quer * w_quer) - w_quer)) / s_sq
	*mata: 		beta = ((s_sq + (w_quer * w_quer) - w_quer) * (w_quer - 1)) / s_sq
		
	mata: 		alpha
	mata: 		beta
	*pause
	
	*mata: 		mata describe
	
	mata: 		st_numscalar("alpha", alpha)
	mata: 		st_numscalar("beta", beta)
	
	restore
	
	**# 	Calculate probabilites for the weights under the beta distribution with alpha and beta
	tempvar 	wi_n
	qui sum		`wgtvar'_norm
	gen 		`wi_n' = `wgtvar'_norm / r(N)
	*pause
	
	tempvar		prob_wgt
	gen 		`prob_wgt' = ibeta(alpha,beta,`wi_n')
	*pause
	
	**#		Compare probabilities to thresholds and set boundaries
	qui sum 	`wgtvar'_norm if `prob_wgt' <= `lowercut' & !missing(`prob_wgt',`wgtvar')
	scalar 		lowerbd = r(max)
	if 			lowerbd == . {
		qui sum 	`wgtvar'_norm
		scalar 		lowerbd = r(min)
	}
	qui sum 	`wgtvar'_norm if `prob_wgt' >= `uppercut' & !missing(`prob_wgt',`wgtvar')
	scalar 		upperbd = r(min)
	if 			upperbd == . {
		qui sum 	`wgtvar'_norm
		scalar 		upperbd = r(max)
	}
	*pause
	
	**# 	Trim weights
	replace 	`wgtvar'_trim = lowerbd if `prob_wgt' <= `lowercut' & !missing(`prob_wgt',`wgtvar')
	replace 	`wgtvar'_trim = upperbd if `prob_wgt' >= `uppercut' & !missing(`prob_wgt',`wgtvar')
	*pause
	
	**# 	Distribute the trimmed weights evenly over the non-trimmed weights
	putmata		orig = `wgtvar'_norm, omit replace
	putmata		trim = `wgtvar'_trim, omit replace
	mata: 		lost = (abs(orig :- trim))' * I
	mata: 		lost
	*pause
	
	mata: 		st_numscalar("lost", lost)
	qui sum 	`wgtvar'_norm if ///
				(`prob_wgt' <= `lowercut' & !missing(`prob_wgt',`wgtvar')) | ///
				(`prob_wgt' >= `uppercut' & !missing(`prob_wgt',`wgtvar'))
	scalar 		lost = lost/r(N)
	replace 	`wgtvar'_trim = `wgtvar'_trim + lost if !( ///
				(`prob_wgt' <= `lowercut' & !missing(`prob_wgt',`wgtvar')) | ///
				(`prob_wgt' >= `uppercut' & !missing(`prob_wgt',`wgtvar')) ///
				)
	
	**#		Normalize trimmed weights again
	qui sum 	`wgtvar'_trim 
	replace 	`wgtvar'_trim = `wgtvar'_trim / r(mean)
	qui sum 	`wgtvar'_trim
	
	**#		Collect information on trimming process
	mat 		trimdok[`round'+1,1] = r(min)
	mat 		trimdok[`round'+1,2] = r(max)
	mat 		trimdok[`round'+1,3] = r(mean)
	mat 		trimdok[`round'+1,4] = r(Var)
}	
	
mat list 		trimdok