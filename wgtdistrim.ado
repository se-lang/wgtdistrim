*! version 0.1.0  10nov2023
program wgtdistrim
    
    version 16.1
    
    syntax varname(numeric) [ if ] [ in ] ///
    , Generate(namelist max=2)            ///
    [                                     ///
        CUToff(real .001)                 ///
        ITERate(integer 10)               ///
        TOLerance(real 0)                 ///
    ]
    
    marksample touse
    quietly replace `touse' = 0 if (`varlist' == 0)
    
    local wgtvar : copy local varlist
    
    typlist_and_varlist_of `generate'
    
    /*
        We now have local macros
        
            wgtvar  <varname> holding untrimmed weights
            typlist  type of new variable from generate(), e.g., float
            varlist  name of new variable from generate()
    */
    
    if ( (`cutoff'<=0) | (`cutoff'>=1) ) ///
        option_invalid cutoff() 125
    
    if (`iterate' < 1) ///
        option_invalid iterate() 125
    
    if (`tolerance' < 0) ///
        option_invalid tolerance() 125
    
    capture assert `wgtvar' >= 0 if `touse' , fast
    if ( _rc ) error 402 // negative weights encountered  
    
    mata : wgtdistrim(       ///
        st_local("wgtvar"),  ///
        st_local("touse"),   ///
        `cutoff',            ///
        `iterate',           ///
        `tolerance',         ///
        st_local("typlist"), ///
        st_local("varlist")  ///
        )
    
end


program typlist_and_varlist_of
    
    syntax newvarname(numeric)
    
    if ("`varlist'" == "") ///
        option_invalid generate() 102
    
    c_local typlist : copy local typlist
    c_local varlist : copy local varlist
    
end


program option_invalid
    
    args option rc
    
    display as err "option `option' invalid"
    exit `rc'
    
end


/*  _________________________________________________________________________
                                                                     Mata  */

version 16.1


mata :


mata set matastrict   on
mata set mataoptimize on


void wgtdistrim(
    
    string scalar wgtvar,
    string scalar touse,
    real   scalar cutoff,
    real   scalar iter,
    real   scalar tolerance,
    string scalar typlist,
    string scalar varlist
    
    )
{
    real colvector w_kt
    real scalar    n
    real scalar    i
    
    
    w_kt = st_data(., wgtvar, touse)
    n    = rows(w_kt)
    
    confirm_obs_and_weights(n, w_kt)
    
    summarize_iteration(0, minmax(w_kt), .)
    
    for (i=1; i<=iter; i++) {
    	
        if (mreldif_w_kt(w_kt,n,cutoff,i) <= tolerance)
            break
        
    }
    
    st_store(., st_addvar(typlist,varlist), touse, w_kt)
}


real scalar mreldif_w_kt(
    
    real colvector w_kt,
    real scalar    n,
    real scalar    cutoff,
    real scalar    iteration
    
    )
{
    real scalar    w_bar
    real scalar    s2
    real scalar    alpha
    real scalar    beta
    real rowvector w_op
    real matrix    K
    real scalar    mreldif
    
    
    w_bar = mean(w_kt)
    s2    = quadcolsum((w_kt:-w_bar):^2/n)
    /*
        (6) in Potter (1990, 227) uses the population variance
        dividing by n, not (n-1).
        
        Perhaps, we should rather use the standard sample variance
        dividing by (n-1), implemented in
        
    s2    = quadvariance(w_i)
    */
    
    alpha = (w_bar*(n*w_bar-1) / (n*s2)) + 2
    beta  = (n*w_bar-1)*(alpha-1)
    
    w_op = 1:/(n*invibetatail(alpha,beta,(cutoff,1-cutoff)))
    
    /*
        The notation below follows Chen et al. (2017, 232)
    */
    
    K = (w_kt:<=w_op[1], w_kt:>=w_op[2])
    
    mreldif = any(K) ? trim_weights(w_kt, w_op, K) : 0
    
    summarize_iteration(iteration, minmax(w_kt), mreldif)
    
    return(mreldif)
}


real scalar trim_weights(
    
    real colvector w_kt,
    real rowvector w_op,
    real matrix    K
    
    )
{
    real colvector w_was
    real scalar    gamma
    
    
    w_was = w_kt
    
    gamma = (
        (quadcolsum(w_kt)-quadsum(K:*w_op)) 
        / 
        quadcolsum((1:-rowsum(K)):*w_kt)
        )
    
    w_kt = gamma:*w_kt
    
    if ( any(K[,1]) )
        w_kt[selectindex(K[,1])] = J(colsum(K[,1]),1,w_op[1])
    
    if ( any(K[,2]) )
        w_kt[selectindex(K[,2])] = J(colsum(K[,2]),1,w_op[2])
    
    return( mreldif(w_was,w_kt) )
}


void summarize_iteration(
    
    real scalar    iteration,
    real rowvector min_max,
    real scalar    mreldif
    
    )
{
	printf("{txt}Iteration %f:", iteration)
    printf("{col 20}{txt}min  = {res}%9.0g", min_max[1])
    printf("{col 40}{txt}max  = {res}%9.0g", min_max[2])
    printf("{col 60}{txt}diff = {res}%9.0g", mreldif)
    printf("\n")
}


void confirm_obs_and_weights(
    
    real scalar    n, 
    real colvector w_i
    
    )
{
    if (n < 2) 
        exit(error(2000+n)) // insufficient observations
    
    /*
        Potter (1990, 227)
        
            f(w) = n (1/nw)^(a+1)*(1 - 1/nw)^(b-1) / B(a,b)
            
                for 1/n <= w <= .
    */
    
    if ( !all((1/n):<=w_i) ) {
        
        errprintf("weights must be greater than %f\n", 1/n)
        exit(459)
        
    }
}


end


exit