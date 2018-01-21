capture program drop FanLiTest1999

program FanLiTest1999, rclass
version 13
//Example of running: FanLiTest1999 eps x1 x2
// Here eps - residuals, x1 and x2 are explanatory variables
// The number of variables can be any

syntax varlist(min=2) 


scalar Lag=1
local Start= 1+Lag


scalar sm = 0 // running sum of what becomes I_n, i.e. statistics
scalar sigma2_hat = 0 // estimate of sigma2
scalar I_n=0 // statistics of interest
scalar Test_stat=0 // test statistic 
scalar pnvar=0 // number of explanatory variables, i.e. regressors



scalar nvar=1 // to be counted number of variables- including one column of errors and several regressors
scalar sd = 1 // standard deviation, start of the product 
local k = 2 // first program argument variable
	
				while "``k''" != "" {
									quietly summarize(``k''), noempty
									//note, that it assumes missing values  '.' as value when computes st.dev.
									scalar sd = sd * r(sd)
									local ++k
									scalar nvar=nvar+1

								}
								
							
								
// compute number of non-missing variables, i.e. excluding "."
// take the first non-error regressor
egen v1_tmp = rownonmiss(`2')
gen cum_amount_tmp = sum(v1_tmp)
scalar N1 = cum_amount_tmp[_N]
drop v1_tmp
drop cum_amount_tmp


scalar pnvar = nvar - 1 //number of explonatary variables, i.e. withouts error terms
scalar hsmooth = ((sd)^(1/(nvar - 1)))*(N1)^(-1/(4 + pnvar)) //bandwith


forval i=`Start'/`=_N' {
// in most cases in stata Start=2 because for lag variables first row is ".", i.e. missing number

		forval j=`Start'/`=_N' {

				if `i' != `j' {
	
				local k = 2 // first program argument variable
				
				scalar product = 1 // product of one-dimensional Kernels
				
				while "``k''" != "" {
	
						MyKernel ``k''[_n+`i'-1] ``k''[_n+`j'-1] hsmooth
						scalar product = product * r(vKernel)						

						local ++k
							}

					scalar  tmp = `1'[_n+`i'-1]*`1'[_n+`j'-1]*product 
					scalar  sm = sm + tmp
					scalar  sigma2_hat=sigma2_hat+tmp*tmp
				
				
				}
				
		}
		

}


 scalar I_n = sm/(N1*N1*hsmooth^pnvar)
 scalar  sigma2_hat = sigma2_hat/(N1*N1*hsmooth^pnvar)
 scalar Test_stat = (N1*hsmooth^(pnvar/2))*I_n/(sigma2_hat^0.5)
 
 //scalar pv_Test_stat= normal(-abs(Test_stat))
  scalar pv_Test_stat= normal(-Test_stat)



return scalar N = N1
return scalar sm = sm
return scalar sd = sd
return scalar nvar	=	nvar
return scalar bandwidth	=	hsmooth
return scalar sigma2_hat = sigma2_hat
return scalar I_n=I_n
return scalar Test_stat = Test_stat
return scalar pv_Test_stat = pv_Test_stat


display as text "Fan and Li (1999) model specification test "
display as text "Under H0 of correct model specification, test statistic is asymptotically N(0,1) "
display as text "Under HA, test statistic diverges to +inf"
display as text "Thus, note that this is a one-sided test"
display as text "Also, note that this test might be undersized in small samples"
display as text "Number of observations:" as result N1 
display as text "Bandwidth (based on Hsiao and Li [2001] rule of thumb): = " as result hsmooth
display as text "Employed kernel: standard normal multivariate probability density function (as in Hsiao and Li [2001])"

display as text "Length of variables(omitting missing elements) = " as result N1
display as text " "
display as text "Test_stat = 	" as result Test_stat
display as text "pv_Test_stat = " as result pv_Test_stat


end

///////////////////////////////////////////////
capture program drop MyKernel
program MyKernel, rclass
// compute one-dimensional Kernel (normal p.d.f.)
	version 13
	args x y hsmooth


	scalar mysum = 0
	scalar vKernel = 0

	scalar vKernel = normalden( (`x' - `y')/`hsmooth')
	return scalar vKernel = vKernel

end

////////////////////////////////////////////////





