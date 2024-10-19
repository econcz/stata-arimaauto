*! version 1.0.7  07oct2024

version 15.1
clear all

loc RS        real scalar
loc RM        real matrix
loc SS        string scalar
loc TM        transmorphic matrix
loc CC        class
loc VV        void

loc version   "version 15.1:"

mata:
mata set matastrict on

`CC' ARIMAAuto extends AssociativeArray                            /* class   */
{
	public:
		`VV' new(), put(), start()
	protected:
		`RS' f_s, f_i, f_sw, f_t, level
		`RM' L, T, MS
		`SS' varlist, ifin, iw, o_hegy, o_dfgls, o_kpss, o_arima, mode, ic
	virtual `RM' get_cv_seas()
}

`VV' ARIMAAuto::new()                                              /* public  */
{
	/*
		                                                                        
		(void) new()                                                            
		                                                                        
		class constructor: default order and values of member variables         
		L, T (optional) and MS should be declared last!                         
	*/
	`TM' tmp, rc
	rc = st_tempname()

	// produce error 1 on any _rc caused by pressing the break key              
	(void) _stata("`version' loc " + rc + " = _rc", 1)
	if (strtoreal(st_local(rc)))                                exit(error(  1))

	// default values of member variables                                       
	/* default IC, lag selection MODE  and confidence LEVEL                   */
	this.level   = /* restricted by critical values in HEGY test   */ c("level")
	this.mode    = /* lag selection for HEGY, DFGLS and KPSS tests */ "maic"
	this.ic      = /* best model selection for ARIMA               */ "aic"
	/* default model and test OPTIONS                                         */
	this.o_hegy  = /* ignored: Mode(), MAXLag(); new option: noGls */ ""
	this.o_dfgls = /* ignored: MAXLag()                            */ ""
	this.o_kpss  = /* ignored: MAXLag()                            */ ""
	this.o_arima = /* ignored: [M][M][A]R[IMA](), ITERate()        */ ""
	/* default FLAGS                                                          */
	this.f_s     = /* consider seasonal models                     */ 1
	this.f_i     = /* consider non-stationary models               */ 1
	this.f_sw    = /* stepwise traversing of the model space       */ 1
	this.f_t     = /* display output from hegy, kpss and arima     */ 0
	/* default HK algorithm LIMITS and MAXIMUM VALUES (1 x 8)                 */
	this.L       = (                                                            
		/* limit of #p in arima(#p,#d,#q)                          */ 5,        
		/* limit of #q in arima(#p,#d,#q)                          */ 5,        
		/* limit of #P in sarima(#P,#D,#Q,#s)                      */ 2,        
		/* limit of #Q in sarima(#P,#D,#Q,#s)                      */ 2,        
		/* limit of inverse characteristic roots                   */ 1/1.001,  
		/* maximum lag order to be used in unit root tests         */ .,        
		/* maximum number of models (rows) in the MODEL SPACE      */ .,        
		/* maximum number of iterations                            */ 100       
	)
	/* default statistics of the unit root TESTS (. x 6)                      */
	this. T      =                                                    J(0,6,.)
	/* default HK algorithm base (4 x 4) for the MODEL SPACE (. x 12)         */
	this.MS      = (                                                            
		/* ARIMA(2,d,2) and ARIMA(2,d,2)(1,D,1)                    */ 2,2,1,1\  
		/* ARIMA(0,d,0) and ARIMA(0,d,0)(0,D,0)                    */ 0,0,0,0\  
		/* ARIMA(1,d,0) and ARIMA(1,d,0)(1,D,0)                    */ 1,0,1,0\  
		/* ARIMA(0,d,1) and ARIMA(0,d,1)(0,D,1)                    */ 0,1,0,1   
	)
}

`VV' ARIMAAuto::start()                                            /* public  */
{
	/*
		                                                                        
		(void) start()                                                          
		                                                                        
		                                                                        
		returns: N/A, updates: MS in the associative array              (. x 12)
	*/
	`RM' ms, c, i, j, k
	`TM' tmp, rc
	rc = st_tempname()

	// general configuration                                                    
	if (exists("varlist")) varlist = get("varlist")                /* varlist */
	else {
		errprintf("Error: results of arima not found\n")
		exit(301)
	}
	if (cols(MS) != 12) {                                          /* models  */
		errprintf("Error: model space undefined\n")
		exit(3001)
	}
	ifin    = exists("ifin") ? (regexm(get("ifin"), "if") ? regexr(get("ifin"), 
	          "if", "if \`touse'") : "if \`touse' "+get("ifin")) : "if \`touse'"
	iw      = exists("iw")   ? get("iw")   : ""
	c       = (ic:==("llf", "aic","sic")) * (10::12)

	// recursive estimation                                                     
	ms = MS                                          /* current model space   */
	MS = J(0,cols(MS),.)                             /* final model space     */
	k  = 0
	if ((tmp=_stata("`version' cap marksample touse", 1)))      exit(error(tmp))
	do {
		/* ARIMA for all IC == . where inverse characteristic roots ≤ L[#]    */
		for(j = 1; j <= rows(ms); j++) {
			// // produce error 1 on any _rc caused by pressing the break key   
			(void) _stata("`version' loc " + rc + " = _rc", 1)
			if (strtoreal(st_local(rc)))                        exit(error(  1))
			// estimate a [S]ARIMA[X] model                                     
			if (_stata("`version' arima " + varlist + " " + ifin + " " + iw    +
			       ", arima("+invtokens(strofreal(ms[j,1..3]))+") "            +
			(f_s ? " sarima("+invtokens(strofreal(ms[j,4..7]))+") "             
			     : "") + (! ms[j,8] ? "noconst " : "")                         +
			       "iter("+strofreal(L[8])+") " + o_arima, 1-(f_t>1))) continue
			if (! st_numscalar("e(converged)"))                        continue
			// store its inverse characteristic roots                           
			if (_stata("`version' cap estat aroots, nogra", 1))        continue
			if (sum(st_matrix("r(Modulus_ar)") :> L[5])                       +
			    sum(st_matrix("r(Modulus_ma)") :> L[5]))               continue
			ms[j,9]      = max((max(st_matrix("r(Modulus_ar)")),                
				                max(st_matrix("r(Modulus_ma)"))))
			// store its LLF and IC                                             
			if (_stata("`version' cap estat ic", 1))                   continue
			if (anyof(st_matrix("r(S)")[(3,5,6)]:==., 1))              continue
			ms[j,10..12] = st_matrix("r(S)")[(3,5,6)]
			// store the model itself                                           
			MS = MS\ms[j,]
			// print message                                                    
			if (! f_t) printf("Model" + strofreal(++k) + ": ARIMA"             +
			       "(" + invtokens(strofreal(ms[j,1..3]), ",") + ")"           +
			(f_s ? "(" + invtokens(strofreal(ms[j,4..7]), ",") + ")" : "")     +
			(ms[j,8] ? " with constant" : " without constant") + "\n")
		}
		if (! rows(MS)) {                            /* if all models failed  */
			errprintf("Error: convergence failure or at least one eigenvalue " +
			          "is at least " + strofreal(L[5]) + "\n")
			exit(430)
		}
		/* stepwise traversing of the MS                                      */
		i = sort(MS, c)[1,]                          /* min IC: MS[1,.]       */
		if (f_sw) {
			ms = abs((J(1,7,0),-1,J(1,cols(i)-8,.)) + i)\(                      
				(((1,0,0,0,0,0\-1,0,0,0,0,0\0,0,1,0,0,0\0,0,-1,0,0,0)\          
				  (1,0,1,0,0,0\-1,0,-1,0,0,0)\                                  
				  (f_s ? (0,0,0,1,0,0\0,0,0,-1,0,0\0,0,0,0,0,1\0,0,0,0,0,-1)\   
				         (0,0,0,1,0,1\0,0,0,-1,0,-1) : J(0,6,.))                
				 ),J(f_s?12:6,6,0)):+(J(f_s?12:6,1,1)#(i[1..8],J(1,4,.))))
			for(j = 3; j <= rows(ms); j++) {         /* limit #p-#Q: <0,L[#]> */
				ms[j,1] = (tmp=ms[j,1]) < 0 ? 0 : (tmp > L[1] ? L[1] : tmp)
				ms[j,3] = (tmp=ms[j,3]) < 0 ? 0 : (tmp > L[2] ? L[2] : tmp)
				ms[j,4] = (tmp=ms[j,4]) < 0 ? 0 : (tmp > L[3] ? L[3] : tmp)
				ms[j,6] = (tmp=ms[j,6]) < 0 ? 0 : (tmp > L[4] ? L[4] : tmp)
			}
			ms = uniqrows(ms)                        /* drop duplicate models */
		}
	} while(f_sw & i[,c] < min(select(MS, MS[,c]:!=i[,c])))
	// final estimation                                                         
	ms = sort(MS, c)[1,]                             /* min IC: MS[1,.]       */
	if ((tmp=_stata("`version' arima " + varlist + " " + ifin + " " + iw       +
	                ", arima("+invtokens(strofreal(ms[1..3]))+") "             +
	         (f_s ? " sarima("+invtokens(strofreal(ms[4..7]))+") " : "")       +
	     (! ms[8] ? " noconst "                                    : "")       +
	                "iter(" + strofreal(L[8]) + ") " + o_arima, 1)))         ///
	                                                            exit(error(tmp))

	// update MS in the associative array                                       
	super.put("MS", MS)
}

`VV' ARIMAAuto::put(`SS' key, `TM' val)                            /* public  */
{
	/*
		                                                                        
		(void) put(key, val)                                                    
		           ---  ---                                                     
		key: *                                                    string (1 x 1)
		val: depvar, ifin, iw, o_*, mode, ic                      string (1 x 1)
		     f_*, level                                           real   (1 x 1)
		     L                                                    string (1 x 5)
		     T (optional)                                         real   (. x 4)
		     MS                                                   string (1 x 2)
		                                                                        
		returns: N/A, updates: * in the associative array                (. x .)
	*/
	`RS' s, D, d, c
	`TM' tmp, rc
	rc = st_tempname()

	// general configuration                                                    
	level = exists("level") & (tmp=get("level")) != . ? tmp : this.level
	if (! anyof(level:==(90,95,99), 1)) {
		errprintf("Error: level must be set to 90, 95 or 99\n")
		exit(198)
	}
	super.put("level", level)
	mode  = exists("mode") & (tmp=get("mode"))   != "" ? tmp : this.mode
	if (! anyof(mode:==("maic","bic","seq"), 1)) {
		errprintf("Error: mode must be maic, bic or seq\n")
		exit(198)
	}
	super.put("mode", mode)
	ic    = exists("ic") & (tmp=get("ic"))       != "" ? tmp : this.ic
	if (! anyof(ic:==("llf","aic","sic"), 1)) {
		errprintf("Error: ic must be llf, aic or sic\n")
		exit(198)
	}
	super.put("ic", ic)

	// all variables except L, T and MS                                         
	if (! anyof(key:==("MS","L"), 1)) super.put(key, val)

	// L                                                                        
	if (key == "L") {
		/* specific configuration                                             */
		if (((tmp=length(tokens(val[1]))) & tmp != 2)                          |
		    ((tmp=length(tokens(val[2]))) & tmp != 2)) {                        
			errprintf("Error: " + (length(tokens(val[1])) ? "max() " : "")     +
			                     (length(tokens(val[2])) ? "mmax() " : "")     +
			          "invalid -- invalid numlist has too few elements\n")
			exit(122)
		}
		if (length(tmp=tokens(val[1]))) L[1..2] = strtoreal(tmp)
		if (length(tmp=tokens(val[2]))) L[3..4] = strtoreal(tmp)
		if (length(val) == 6          ) L[5..8] = strtoreal(val[3..6])
		/* update L in the associative array                                  */
		super.put("L",  L )
	}

	// T and MS                                                                 
	if (key == "MS") {
		/* specific configuration                                             */
		if (((tmp=length(tokens(val[1]))) & tmp != 3)                          |
		    ((tmp=length(tokens(val[2]))) & tmp != 4)) {                        
			errprintf("Error: " + (length(tokens(val[1])) ? "arima() " : "")   +
			                     (length(tokens(val[2])) ? "sarima() " : "")   +
			          "invalid -- invalid numlist has too few elements\n")
			exit(122)
		}
		if (exists("varlist")) varlist = tokens(get("varlist"))[1] /* depvar  */
		else {
			errprintf("Error: results of arima not found\n")
			exit(301)
		}
		ifin    = exists("ifin")    ? get("ifin")               : ""
		/* AssociativeArray -> member variables                               */
		o_hegy  = exists("o_hegy")  ? regexr(regexr(regexr(get("o_hegy"),       
		                              "maxl[ag]*\([^(]+\)", ""),                
		                              "^[no]*g[ls]*| [no]*g[ls]*", ""),         
		                              "m[deo]*\([^(]+\)", "")   : this.o_hegy
		                              /* ignored: Gls, Mode(), MAXLag()       */
		o_dfgls = exists("o_dfgls") ? regexr(get("o_dfgls"),                    
		                              "maxl[ag]*\([^(]+\)", "") : this.o_dfgls
		                              /* ignored: MAXLag()                    */
		o_kpss  = exists("o_kpss")  ? regexr(get("o_kpss"),                     
		                              "maxl[ag]*\([^(]+\)", "") : this.o_kpss
		                              /* ignored: MAXLag()                    */
		o_arima = exists("o_arima") ? regexr(regexr(regexr(regexr(regexr(regexr(
		                              regexr(get("o_arima"), "arima\([^(]+\)",  
		                              ""), "sarima\([^(]+\)", ""),              
		                              "ar\([^(]+\)", ""), "ma\([^(]+\)", ""),   
		                              "mar\([^(]+\)", ""), "mma\([^(]+\)", ""), 
		                              "iterate\([^(]+\)", "") : this.o_arima
		                              /* ignored: [MMA]R[IMA](), ITERate()    */
		f_s     = exists("f_s")     ? get("f_s")                : this.f_s
		f_i     = exists("f_i")     ? get("f_i")                : this.f_i
		f_sw    = exists("f_sw")    ? get("f_sw")               : this.f_sw
		f_t     = exists("f_t")     ? get("f_t")                : this.f_t
		L       = exists("L" )      ? get("L" )                 : this.L
		T       = exists("T" )      ? get("T" )                 : this.T
		MS      = exists("MS")      ? get("MS")                 : this.MS
		/* seasonal period and difference, #s and #D in sarima(#P,#D,#Q,#s)   */
		s = D = 0
		if (f_s & ! length(tokens(val[2]))) {
			// estimate #s, disable seasonal models for yearly data             
			if ((tmp=_stata("`version' cap tsset, noq", 1)))    exit(error(tmp))
			f_s = s = (st_global("r(unit1)"):==("q","m")) * (4\12)
			super.put("f_s", f_s != 0)
			// estimate #D if #s > 0                                            
			if (f_i & s) {
				if ((tmp=_stata("`version' hegy " + varlist + " " + ifin       +
				    ", " + (! regexm(o_hegy, "^nog| nog") ? "g " : " ")        +
				    "m("+mode+") "                                             +
				    (L[6] < . ? "maxl("+strofreal(L[6])+") " : "")             +
				    o_hegy, 1-(f_t>0))))                        exit(error(tmp))
				T = T\((                                                        
					.,                               /* seasonal unit root    */
					.,                               /* number of lags        */
					st_numscalar("r(t2)")),                                     
					get_cv_seas(s, (! regexm(o_hegy, "^nog| nog") ? "g " : "") +
					            o_hegy)[2,])
				if (T[1,3] > T[1,cols(T)+trunc((100-level)/5-2)]) D = 1
				/* H0: seasonal unit root, t_S/2 > t_cv                       */
			} else if (f_i & st_global("r(unit1)") != ".") {
				errprintf(                                                      
					"hegy must be used with monthly or quarterly data\n"       +
					"please re-tsset the data or add a sarima() option\n"       
				)
				exit(198)
			}
		} else if (f_s) D = strtoreal(tokens(val[2]))[2]
		/* trend difference, #d in arima(#p,#d,#q)                            */
		d = 0
		if (f_i & ! length(tokens(val[1]))) {
			mode = (tmp=invtokens((mode:==("bic","seq")):*("sic","opt"))=="")  ?
				   tmp : mode
			d--
			do {
				// produce error 1 on any _rc caused by pressing the break key  
				(void) _stata("`version' loc " + rc + " = _rc", 1)
				if (strtoreal(st_local(rc)))                    exit(error(  1))
				// estimate #d recursively                                      
				d++
				if ((tmp=_stata("`version' dfgls "                             +
				    "D" + strofreal(d) + ".S" + strofreal(s) + "." + varlist   +
				    " " + ifin + ", "                                          +
				    (L[6] < . ? "maxl("+strofreal(L[6])+") " : "")             +
				    o_dfgls, 1-(f_t>0))))                       exit(error(tmp))
				T = T\((                                                        
					.,                               /* unit root             */
					st_numscalar("r("+mode+"lag)")   /* number of lags        */
				),(tmp=(st_matrix("r(cvalues)") == J(0,0,.)                     
					? J(st_numscalar("r("+mode+"lag)")+1,5,.)                   
					: st_matrix("r(cvalues)")                                   
				))[                                                             
					rows(tmp)-st_numscalar("r("+mode+"lag)")+1,                 
					2..5                                                        
				])
				if ((tmp=_stata("`version' kpss "                              +
				    "D" + strofreal(d) + ".S" + strofreal(s) + "." + varlist   +
				    " " + ifin + ", "                                          +
				    "maxl("+strofreal(st_numscalar("r("+mode+"lag)"))+") "     +
				    o_kpss, 1-(f_t>0))))                        exit(error(tmp))
				T = T\((                                                        
					.,                               /* unit root             */
					T[rows(T),2],                    /* number of lags        */
					st_numscalar("r(kpss" + strofreal(T[rows(T),2]) + ")")      
				),(! regexm(o_kpss, "^not| not") ? (0.216,0.146,0.119) :        
					                               (0.739,0.463,0.347)))
			} while(T[rows(T)-1,3]>T[rows(T)-1,cols(T)+trunc((100-level)/5-2)] &
			        /* H0: unit root,    tau > t_cv                          */ 
			        T[rows(T),  3]>T[rows(T),  cols(T)+trunc((100-level)/5-2)]  
			        /* H0: stationarity, LM  < LM_cv                         */)
		} else if (f_i) d = strtoreal(tokens(val[1]))[2]
		/* lags, statistics and critical values of unit root tests, T         */
		T[,1] = T[,3] :> T[,cols(T)+trunc((100-level)/5-2)]
		// update T in the associative array                                    
		super.put("T", T)
		/* initial model space with|without constant, MS                      */
		c   = (length((tokens(val[1]),tokens(val[2]))) == 7                     
				? (! regexm(o_arima, "^noc| noc") ? 1 : 0)                      
				: (d + D < 2 ? 1 : 0)                                           
		)
		// stepwise traversing of the MS                                        
		tmp = (f_s ? (MS,J(rows(MS),1,d),J(rows(MS),1,D),J(rows(MS),1,s)        
			         )[,(1,5,2,3,6,4,7)]                                        
			       : (MS[,1..2],J(rows(MS),1,d))[,(1,3,2)],                     
			         J(rows(MS),4,0)                                            
		)
		MS  = (! length(tokens(invtokens(val)))      /* default model space   */
			? tmp                                                               
			: (length(tokens(invtokens(val))) <= 4   /* user-defined model(s) */
				? (! length(tokens(val[1]))                                     
					? tmp[,1..3], (J(rows(tmp),1,1)#strtoreal(tokens(val[2])))  
					: (J(rows(tmp),1,1)#strtoreal(tokens(val[1]))),tmp[,4..7]   
				  )                                                             
				: strtoreal((tokens(val[1]),tokens(val[2])))                    
			  )                                                                 
		)
		// bulk estimation of the MS                                            
		if (! f_sw) {
			MS  = (                                                             
				(! length(tokens(invtokens(val)))                               
					? J(0,7,.)                                                  
					: MS                                                        
				)\                                   /* user-defined model(s) */
				(f_s ? ((tmp=(((0::L[1])#J(L[2]+1,1,1),                         
					     J(L[1]+1,1,1)#(0::L[2]))#J(L[3]+1,1,1),                
					     J((L[1]+1)*(L[2]+1),1,1)#(0::L[3]))#J(L[4]+1,1,1),     
					     J((L[1]+1)*(L[2]+1)*(L[3]+1),1,1)#(0::L[4])),          
					    J(rows(tmp),1,d),J(rows(tmp),1,D),J(rows(tmp),1,s)      
					   )[,(1,5,2,3,6,4,7)]                                      
					 : ((tmp=(0::L[1])#J(L[2]+1,1,1),J(L[1]+1,1,1)#(0::L[2])),  
					    J(rows(tmp),1,d))[,(1,3,2)],J(rows(tmp),4,0)            
				)                                                               
			)[|1,.\L[7],.|]                          /* limit MS if specified */
		}
		MS  = MS,J(rows(MS),1,c),J(rows(MS),4,.)
		// update MS in the associative array                                   
		super.put("MS", MS)
	}
}

`RM' ARIMAAuto::get_cv_seas(`RS' s, `SS' o)                        /* virtual */
{
	/*
		--                                                                      
		cv = `RM' get_cv_seas(s, o)                                             
		                      -  -                                              
		HEGY critical values at 1%, 5% and 10% significance levels, calculated  
		from the response surface equation:                                     
		q^p(T) = θ_∞^p + θ_1^p * T^−1 + θ_2^p * T^−2 + θ_3^p * T^−3 + ε         
		with coefficients provided by Barrio Castro, Bodnar and Sansó (2015)    
		                                                                        
		returns: (t_0\t_S/2\F_K\F_SEAS\F_ALL x 0.01,0.05,0.10)           (5 x 3)
	*/
	`RM' N, theta
	`TM' tmp
	N   = st_numscalar("e(N)")                                        /* _N   */
	tmp = (! regexm(o, "^g| g") ? "ols" : "gls") + strofreal(s) + " "          +
	      (regexm(o, "det\([^(]+\)") ? regexs(0) : "det(seas)")       /* CASE */

	// response surface equation coefficients                                   
	/*                   θ_∞^p        θ_1^p        θ_2^p        θ_3^p         */
	/* None   */ if (tmp == "ols4 det(none)") theta = (
		/* t_0    0.01 */  -2.5676678,  0.47516783,  -2.0866101,   11.057007\   
		/*        0.05 */  -1.9411654,  0.65898524,  -1.3763645,   11.727207\   
		/*        0.10 */  -1.6168177,  0.65676244, -0.67963198,   6.1568001\   
		/* t_S/2  0.01 */  -2.5677929,  0.39461038, -0.57613774,   3.9353396\   
		/*        0.05 */  -1.9410209,  0.56672883,   1.1513010,  -5.4928701\   
		/*        0.10 */  -1.6167719,  0.59752691,  0.68462318,  -1.8546284\   
		/* F_K    0.01 */   4.7298820,   1.2719035,  -1.2462890,   15.707069\   
		/*        0.05 */   3.1105441, -0.47461635, -0.58267156,   7.0978537\   
		/*        0.10 */   2.4094103, -0.88414709,   1.7175107,  -9.3520629\   
		/* F_SEAS 0.01 */   3.9360076,   2.1138295,   1.5101301,   2.2413873\   
		/*        0.05 */   2.7444305,  0.30960883, -0.99512150,   7.6669009\   
		/*        0.10 */   2.2154677, -0.26964074,  0.67703277,  -2.6692922\   
		/* F_ALL  0.01 */   3.4810704,   2.8444620, -0.73302974,   23.059449\   
		/*        0.05 */   2.5212908,  0.78682142, -0.17383853,   3.6302940\   
		/*        0.10 */   2.0866174,  0.14362515,  0.66683400,  -3.1813570    
	)
	/* CASE 1 */ if (tmp == "ols4 det(const)") theta = (
		/* t_0    0.01 */  -3.4279930, -0.58122008,   2.5521031,  -28.067062\   
		/*        0.05 */  -2.8602236,  0.23170727,   1.0513884,  -9.7978441\   
		/*        0.10 */  -2.5660487,  0.49846225,  0.50061794,  -4.2384347\   
		/* t_S/2  0.01 */  -2.5654843,  0.33278964,   2.0609683,  -11.813474\   
		/*        0.05 */  -1.9409957,  0.64407136,  0.19919366,   1.7642012\   
		/*        0.10 */  -1.6169878,  0.66476432,-0.058310126,   3.1426989\   
		/* F_K    0.01 */   4.7324106, -0.10668131,   1.2958689,   28.570727\   
		/*        0.05 */   3.1101654,  -1.3008103,-0.078482600,   18.377518\   
		/*        0.10 */   2.4091106,  -1.4663395,  0.21079267,   13.804404\   
		/* F_SEAS 0.01 */   3.9340198,   1.3414927,  -1.2242778,   37.889254\   
		/*        0.05 */   2.7441775, -0.31384041, -0.85485001,   16.809288\   
		/*        0.10 */   2.2145428, -0.68409937, -0.85121981,   13.143561\   
		/* F_ALL  0.01 */   4.3786638,   4.8517259,  -6.4152528,   81.284739\   
		/*        0.05 */   3.3069139,   1.6576716,  -3.9239019,   37.871160\   
		/*        0.10 */   2.8079755,  0.60790886,  -1.7231322,   15.977384    
	)
	/* CASE 2 */ if (tmp == "ols4 det(trend)") theta = (
		/* t_0    0.01 */  -3.9578302, -0.95063027, 0.045280054,  -18.612543\   
		/*        0.05 */  -3.4096333, 0.057242292,  0.62360541,  -9.3269809\   
		/*        0.10 */  -3.1271451,  0.48224946, -0.52447487,  0.79311180\   
		/* t_S/2  0.01 */  -2.5667639,  0.17693653,  0.69950518,  -5.3095400\   
		/*        0.05 */  -1.9411076,  0.44326384,  0.85436876,  -3.2174850\   
		/*        0.10 */  -1.6167623,  0.46959333,  0.64073955,  -1.1257551\   
		/* F_K    0.01 */   4.7315810,  -1.5950089,   6.5630488,   32.858942\   
		/*        0.05 */   3.1104266,  -2.4251299,   6.5757176,  -2.6263650\   
		/*        0.10 */   2.4090477,  -2.3224461,   5.2016427,  -3.3865231\   
		/* F_SEAS 0.01 */   3.9358326,  0.47449633,  0.54970285,   59.866301\   
		/*        0.05 */   2.7450620,  -1.0033942,   3.0899005,   9.5584463\   
		/*        0.10 */   2.2156122,  -1.2923094,   3.4962851,  -1.1475224\   
		/* F_ALL  0.01 */   5.2504668,   6.3605264,  0.54897335,   83.532860\   
		/*        0.05 */   4.0936030,   2.4357024,  -1.1958470,   36.041727\   
		/*        0.10 */   3.5482854,   1.0664806, -0.60926013,   17.885635    
	)
	/* CASE 3 */ if (tmp == "ols4 det(seas)") theta = (
		/* t_0    0.01 */  -3.4297763,  0.36001376, 0.048111637,  -26.143721\   
		/*        0.05 */  -2.8606730,  0.97223479,  0.26045241,  -9.2990100\   
		/*        0.10 */  -2.5665713,   1.1718978, -0.17083371,  -1.9030244\   
		/* t_S/2  0.01 */  -3.4286971,  0.39772945,  -1.8735404,  -12.559216\   
		/*        0.05 */  -2.8616439,   1.0286748, -0.81450059,  -3.3678362\   
		/*        0.10 */  -2.5669902,   1.1774820, -0.18065850,  -2.3249076\   
		/* F_K    0.01 */   8.8019274,   3.3538848,   14.277935,   70.724390\   
		/*        0.05 */   6.6424614, -0.92667897,   2.4713914,   34.417913\   
		/*        0.10 */   5.6266552,  -2.1387320,-0.065003639,   22.338578\   
		/* F_SEAS 0.01 */   7.5396048,   7.5991821,   8.7130426,   104.98413\   
		/*        0.05 */   5.9104902,   1.9392848,   5.3296319,   18.793816\   
		/*        0.10 */   5.1271615,  0.32007337, -0.95729121,   25.507557\   
		/* F_ALL  0.01 */   6.8331686,   10.088421,   11.745679,   108.79019\   
		/*        0.05 */   5.4859552,   4.2840335,   1.1083491,   45.083524\   
		/*        0.10 */   4.8355430,   2.0705299, -0.86649686,   26.985332    
	)
	/* CASE 4 */ if (tmp == "ols4 det(strend)") theta = (
		/* t_0    0.01 */  -3.9591996,0.0092665103,  -3.2784932,  -20.098089\   
		/*        0.05 */  -3.4091581,  0.76224153,  0.75572437,  -14.450850\   
		/*        0.10 */  -3.1264299,   1.1007609,  0.40405671,  -3.7484708\   
		/* t_S/2  0.01 */  -3.4288406,  0.43970705,  -2.9928732,  -7.7696786\   
		/*        0.05 */  -2.8613813,   1.0097287, -0.86423766,  -1.1243402\   
		/*        0.10 */  -2.5670798,   1.1690247, -0.48726385,   1.8289946\   
		/* F_K    0.01 */   8.8079433,   1.3768422,   27.927537,   43.673534\   
		/*        0.05 */   6.6484906,  -2.4078609,   13.518526,  -6.3163026\   
		/*        0.10 */   5.6313541,  -3.1794977,   5.3193418,   5.0928481\   
		/* F_SEAS 0.01 */   7.5430500,   6.4424174,   17.852255,   94.686737\   
		/*        0.05 */   5.9137395,   1.1626669,   11.496638,  -3.2857448\   
		/*        0.10 */   5.1309255, -0.38618742,   4.8934131,   2.3557045\   
		/* F_ALL  0.01 */   7.6191556,   12.071510,   16.260612,   157.81323\   
		/*        0.05 */   6.2116903,   5.0771543,   9.0521913,   24.511814\   
		/*        0.10 */   5.5243556,   2.7015249,   1.4485982,   23.006445    
	)
	/* CASE 5 */ if (tmp == "ols4 det(mult)") theta = (
		/* t_0    0.01 */  -3.9588101, -0.26099430,  -6.8041457,  -28.554100\   
		/*        0.05 */  -3.4103017,  0.60746397,  -2.5274287,  -11.745063\   
		/*        0.10 */  -3.1270704,  0.92161079,  -1.8381578,  0.62060474\   
		/* t_S/2  0.01 */  -3.9592924, -0.20648349,  -6.9819350,  -31.230843\   
		/*        0.05 */  -3.4104150,  0.59059676,  -2.1206633,  -13.207976\   
		/*        0.10 */  -3.1275135,  0.90402932,  -1.2230855,  -2.8750693\   
		/* F_K    0.01 */   12.207936,   9.5227497,   22.455666,   358.76220\   
		/*        0.05 */   9.7445908,   1.5427798,   9.9435242,   98.777146\   
		/*        0.10 */   8.5712092, -0.98585619,   2.4692880,   52.913858\   
		/* F_SEAS 0.01 */   10.750947,   16.028405,   21.898973,   446.20257\   
		/*        0.05 */   8.8662471,   6.7413723,   10.228860,   135.52841\   
		/*        0.10 */   7.9494832,   3.3764254,   3.6260652,   67.391514\   
		/* F_ALL  0.01 */   9.9188917,   20.406505,   17.886062,   501.64223\   
		/*        0.05 */   8.3530514,   10.128553,   9.8759863,   159.10682\   
		/*        0.10 */   7.5807963,   6.4076414,-0.039221509,   101.85158    
	)
	/* CASE 1 */ if (tmp == "gls4 det(const)") theta = (
		/* t_0    0.01 */  -2.6017744,  -11.576394,   103.59152,  -396.96190\   
		/*        0.05 */  -1.9890882,  -13.271645,   124.01966,  -460.71773\   
		/*        0.10 */  -1.6709987,  -14.719768,   140.06090,  -515.29043\   
		/* t_S/2  0.01 */  -2.5633008, -0.14334257,  -2.9538673,   14.146696\   
		/*        0.05 */  -1.9396597,  0.29797905,  -3.7835909,   21.564794\   
		/*        0.10 */  -1.6154096,  0.33501888,  -2.5383989,   13.526894\   
		/* F_K    0.01 */   4.7455076, 0.088391710,   9.2387492,  -4.5591302\   
		/*        0.05 */   3.1157636, -0.91730683,   2.3154692,  0.16499043\   
		/*        0.10 */   2.4101268, -0.98312773, -0.60507862,   11.350452\   
		/* F_SEAS 0.01 */   3.9375752,   2.0481440,   4.6598029,   8.3126654\   
		/*        0.05 */   2.7463841,  0.22557615,   2.1941347,  0.12343295\   
		/*        0.10 */   2.2153354, -0.21665313,   1.2369000, -0.57292537\   
		/* F_ALL  0.01 */   3.4591416,   13.786047,  -56.194826,   183.58524\   
		/*        0.05 */   2.5003904,   10.751192,  -52.532820,   129.29991\   
		/*        0.10 */   2.0644772,   9.6646754,  -51.480407,   116.29608    
	)
	/* CASE 2 */ if (tmp == "gls4 det(trend)") theta = (
		/* t_0    0.01 */  -3.4294600,  -9.7579046,   82.481958,  -310.64282\   
		/*        0.05 */  -2.8745889,  -10.209661,   100.40630,  -379.15980\   
		/*        0.10 */  -2.5895772,  -10.773139,   112.83004,  -431.11607\   
		/* t_S/2  0.01 */  -2.5648603, -0.72207957,  -1.6630053,   5.1855040\   
		/*        0.05 */  -1.9396038, -0.25439051,  -1.2392589,   9.0240309\   
		/*        0.10 */  -1.6157063, -0.11295406,  -1.0298877,   7.4431335\   
		/* F_K    0.01 */   4.7319420,  0.43404756,   8.7469067,   23.626831\   
		/*        0.05 */   3.1124684,  -1.0782544,   5.9701458,   2.9993925\   
		/*        0.10 */   2.4105288,  -1.2820898,   4.6041636, -0.48669346\   
		/* F_SEAS 0.01 */   3.9359790,   2.1304154,   13.614065,  -16.254418\   
		/*        0.05 */   2.7454947,  0.34711481,   4.7715588,   5.4579053\   
		/*        0.10 */   2.2152630, -0.15978846,   3.8089946,  -1.9892152\   
		/* F_ALL  0.01 */   4.4016081,   20.771288,  -121.55971,   503.01758\   
		/*        0.05 */   3.3400212,   17.111022,  -131.87092,   508.37418\   
		/*        0.10 */   2.8472804,   15.794686,  -137.36540,   528.29785    
	)
	/* CASE 3 */ if (tmp == "gls4 det(seas)") theta = (
		/* t_0    0.01 */  -2.5990665,  -12.047061,   85.687973,  -313.47376\   
		/*        0.05 */  -1.9838830,  -13.928391,   115.12810,  -418.81782\   
		/*        0.10 */  -1.6671559,  -15.262401,   132.19935,  -478.62541\   
		/* t_S/2  0.01 */  -2.5974222,  -12.235628,   89.523079,  -333.16689\   
		/*        0.05 */  -1.9841934,  -13.936898,   115.37681,  -419.67485\   
		/*        0.10 */  -1.6677888,  -15.228615,   131.54416,  -474.62520\   
		/* F_K    0.01 */   4.7198851,   26.092445,  -76.364389,   193.84426\   
		/*        0.05 */   3.0977204,   22.073034,  -68.871248,   99.894313\   
		/*        0.10 */   2.3932811,   20.125496,  -59.807989,   47.881951\   
		/* F_SEAS 0.01 */   3.8908411,   32.625569,  -95.725256,   261.93441\   
		/*        0.05 */   2.7031125,   27.816240,  -96.678618,   197.87319\   
		/*        0.10 */   2.1750663,   25.431190,  -90.365416,   151.30750\   
		/* F_ALL  0.01 */   3.4294484,   34.789140,  -100.53054,   300.41199\   
		/*        0.05 */   2.4707113,   29.916552,  -103.76200,   231.20254\   
		/*        0.10 */   2.0364976,   27.593827,  -100.52671,   190.05771    
	)
	/* CASE 4 */ if (tmp == "gls4 det(strend)") theta = (
		/* t_0    0.01 */  -3.4240594,  -10.555798,   66.839623,  -252.60934\   
		/*        0.05 */  -2.8694267,  -10.939314,   87.891982,  -332.34765\   
		/*        0.10 */  -2.5859015,  -11.344668,   99.533099,  -378.36706\   
		/* t_S/2  0.01 */  -2.6002007,  -12.738076,   87.668208,  -326.35103\   
		/*        0.05 */  -1.9843278,  -14.491218,   117.00689,  -429.04365\   
		/*        0.10 */  -1.6677164,  -15.740811,   133.48634,  -483.72754\   
		/* F_K    0.01 */   4.7204106,   26.102737,  -70.902319,   196.40945\   
		/*        0.05 */   3.0925571,   22.430856,  -70.177911,   106.56658\   
		/*        0.10 */   2.3878910,   20.555888,  -63.130510,   60.869493\   
		/* F_SEAS 0.01 */   3.8906697,   33.381043,  -90.309043,   273.66357\   
		/*        0.05 */   2.7004043,   28.523640,  -94.299930,   194.36907\   
		/*        0.10 */   2.1718384,   26.098291,  -89.188768,   145.03976\   
		/* F_ALL  0.01 */   4.3652945,   41.380478,  -148.50333,   602.59639\   
		/*        0.05 */   3.3081417,   35.649826,  -161.41505,   543.91920\   
		/*        0.10 */   2.8154926,   33.249670,  -166.76296,   535.00317    
	)
	/* CASE 5 */ if (tmp == "gls4 det(mult)") theta = (
		/* t_0    0.01 */  -3.4217626,  -13.179493,   71.560370,  -301.72588\   
		/*        0.05 */  -2.8684513,  -13.120393,   91.332760,  -366.27019\   
		/*        0.10 */  -2.5849263,  -13.352406,   103.26190,  -410.92894\   
		/* t_S/2  0.01 */  -3.4250965,  -13.019092,   69.285867,  -291.10874\   
		/*        0.05 */  -2.8703909,  -13.034822,   90.663643,  -366.62196\   
		/*        0.10 */  -2.5865534,  -13.274299,   102.48787,  -409.96937\   
		/* F_K    0.01 */   8.6583190,   54.666891,  -108.73738,   537.09139\   
		/*        0.05 */   6.5914446,   47.231376,  -173.72113,   638.06151\   
		/*        0.10 */   5.6307484,   44.497686,  -199.79752,   691.47889\   
		/* F_SEAS 0.01 */   7.5183670,   64.231225,  -183.85832,   1080.1022\   
		/*        0.05 */   5.9550030,   56.090770,  -246.30303,   1090.4206\   
		/*        0.10 */   5.2153158,   52.660717,  -268.92965,   1100.1364\   
		/* F_ALL  0.01 */   6.8682745,   67.762763,  -215.13830,   1312.4729\   
		/*        0.05 */   5.5787523,   59.621558,  -273.76368,   1277.0575\   
		/*        0.10 */   4.9624761,   56.003649,  -294.43467,   1265.3690    
	)
	/* None   */ if (tmp == "ols12 det(none)") theta = (
		/* t_0    0.01 */  -2.5675401,   1.3094475,  -2.3015686,   9.8562761\   
		/*        0.05 */  -1.9417622,   1.0719636, -0.51224492,  0.58055632\   
		/*        0.10 */  -1.6175985,  0.96263173,  -1.0662905,   5.0677654\   
		/* t_S/2  0.01 */  -2.5664128,   1.2292061,  -1.2555395,   6.2523839\   
		/*        0.05 */  -1.9407732,  0.99387457,  0.34360137,  -1.5542122\   
		/*        0.10 */  -1.6170914,  0.91633645, -0.54376168,   4.0409648\   
		/* F_K    0.01 */   4.7322746,  -2.5629962,  -1.1271589,   13.311456\   
		/*        0.05 */   3.1095037,  -2.0146195,  -1.7868307,   12.504621\   
		/*        0.10 */   2.4070460,  -1.6520329,  -2.0280066,   12.957738\   
		/* F_SEAS 0.01 */   2.3456470,  0.28491852,  -1.0546447,   7.1881190\   
		/*        0.05 */   1.8775094, -0.13830753,  -1.2240437,   5.5672143\   
		/*        0.10 */   1.6550286, -0.29355936,  -1.0243531,   4.1683709\   
		/* F_ALL  0.01 */   2.2890432,  0.40020356, -0.62816805,   5.9506023\   
		/*        0.05 */   1.8481484,-0.080911444, -0.23332386,  0.60589613\   
		/*        0.10 */   1.6369129, -0.21882533, -0.72160961,   2.8142843    
	)
	/* CASE 1 */ if (tmp == "ols12 det(const)") theta = (
		/* t_0    0.01 */  -3.4307657,   1.2029985,  -2.1075018,   4.2503217\   
		/*        0.05 */  -2.8610683,   1.1901521, -0.84641464,   1.6992724\   
		/*        0.10 */  -2.5674250,   1.2274372,  -1.8215069,   7.6229947\   
		/* t_S/2  0.01 */  -2.5658514,   1.1654154,  0.72217674,  -3.4849549\   
		/*        0.05 */  -1.9400205,  0.96591540,   1.2559339,  -5.9754130\   
		/*        0.10 */  -1.6158989,  0.85108610,  0.94671298,  -3.8560344\   
		/* F_K    0.01 */   4.7347231,  -3.1429516,   3.0630371,  -6.3132977\   
		/*        0.05 */   3.1113755,  -2.4513097,   1.5869998,  -2.5064560\   
		/*        0.10 */   2.4088408,  -2.0257630,   1.2244484,  -1.9689176\   
		/* F_SEAS 0.01 */   2.3447725, 0.025158169,  0.98280235,  -3.6657225\   
		/*        0.05 */   1.8780538, -0.34713990, -0.21384611,   1.8730848\   
		/*        0.10 */   1.6550581, -0.45109894, -0.61094604,   3.5495146\   
		/* F_ALL  0.01 */   2.5342473,  0.57016652, -0.90487963,   10.131331\   
		/*        0.05 */   2.0695665,-0.019938837, -0.62626988,   3.7168939\   
		/*        0.10 */   1.8454335, -0.23289325, -0.42314075,   1.5754821    
	)
	/* CASE 2 */ if (tmp == "ols12 det(trend)") theta = (
		/* t_0    0.01 */  -3.9549846,   1.0122235,   3.4984414,  -30.675131\   
		/*        0.05 */  -3.4083525,   1.2863159,  0.51188055,  -6.9046369\   
		/*        0.10 */  -3.1255590,   1.3222249,  0.37964221,  -4.9865777\   
		/* t_S/2  0.01 */  -2.5662560,   1.1300491,  0.40339422,  -3.0012104\   
		/*        0.05 */  -1.9411086,  0.98045594, -0.15456280,   2.9056916\   
		/*        0.10 */  -1.6165843,  0.84099081, 0.091448416,   1.7285454\   
		/* F_K    0.01 */   4.7319190,  -3.4930734,   1.4612942,   13.170863\   
		/*        0.05 */   3.1100357,  -2.6862001,  0.70787752,   7.7103220\   
		/*        0.10 */   2.4080743,  -2.2397875,   1.1979622,   1.7845587\   
		/* F_SEAS 0.01 */   2.3445369, -0.11785483, -0.63097494,   10.633739\   
		/*        0.05 */   1.8777836, -0.51559365,-0.067879271,   3.9232476\   
		/*        0.10 */   1.6550991, -0.61494772, -0.33241107,   5.1974024\   
		/* F_ALL  0.01 */   2.7887716,  0.98176025,  -5.6804247,   39.999899\   
		/*        0.05 */   2.3078427, 0.057863104,  -1.3829160,   11.691416\   
		/*        0.10 */   2.0734750, -0.22504063, -0.44292362,   4.0147672    
	)
	/* CASE 3 */ if (tmp == "ols12 det(seas)") theta = (
		/* t_0    0.01 */  -3.4305843,   2.3483579,  -3.6773595,   5.3848125\   
		/*        0.05 */  -2.8622944,   2.2297365,  -2.8249245,   8.4082381\   
		/*        0.10 */  -2.5677525,   2.0958616,  -2.4516818,   10.745195\   
		/* t_S/2  0.01 */  -3.4305505,   2.3500403,  -3.6558281,   4.4179047\   
		/*        0.05 */  -2.8606026,   2.0657875,  0.15886793,  -6.0454210\   
		/*        0.10 */  -2.5655316,   1.9159787,  0.87515889,  -6.1630405\   
		/* F_K    0.01 */   8.8059579,  -10.372729,   13.962657,-0.026354168\   
		/*        0.05 */   6.6439349,  -8.9382092,   7.0052185,  0.87993246\   
		/*        0.10 */   5.6291142,  -8.0390517,   5.5872678,  -5.0031789\   
		/* F_SEAS 0.01 */   5.1879385,   1.8400227,  0.83055765,   28.431909\   
		/*        0.05 */   4.4703393,  0.12985135,  -1.4541350,   15.043602\   
		/*        0.10 */   4.1108956, -0.55490497,  -2.3373856,   11.801312\   
		/* F_ALL  0.01 */   5.0832624,   2.4113652,   2.7749309,   15.916600\   
		/*        0.05 */   4.4050600,  0.57786059, -0.43861369,   8.9351938\   
		/*        0.10 */   4.0639292, -0.15341110,  -1.9719314,   10.047546    
	)
	/* CASE 4 */ if (tmp == "ols12 det(strend)") theta = (
		/* t_0    0.01 */  -3.9559888,   2.3216257,  -2.0198937,  -4.5424273\   
		/*        0.05 */  -3.4088741,   2.3214439, -0.85205799,  0.88782690\   
		/*        0.10 */  -3.1255230,   2.2260245,  0.36176910,  -1.8467394\   
		/* t_S/2  0.01 */  -3.4282036,   2.1868232,  -1.0145889,  -6.0207827\   
		/*        0.05 */  -2.8612384,   2.1208895, -0.96760247,   1.7471590\   
		/*        0.10 */  -2.5661124,   1.9492966,  0.19192629,  -1.2086890\   
		/* F_K    0.01 */   8.8084506,  -11.014045,   17.130188,  -7.8454988\   
		/*        0.05 */   6.6445821,  -9.2838275,   7.5795355,  0.43730924\   
		/*        0.10 */   5.6294078,  -8.3078736,   5.8984249,  -5.4548148\   
		/* F_SEAS 0.01 */   5.1872176,   1.6501231,   1.3986009,   22.779383\   
		/*        0.05 */   4.4699058,0.0051320780,  -2.0261104,   15.795071\   
		/*        0.10 */   4.1118809, -0.77106095, -0.95479506,   1.6726229\   
		/* F_ALL  0.01 */   5.3142769,   2.6984230,   2.4181515,   19.104695\   
		/*        0.05 */   4.6223044,  0.79321478,  -2.1911885,   17.976482\   
		/*        0.10 */   4.2759089,-0.099684413,  -1.6781388,   5.9888628    
	)
	/* CASE 5 */ if (tmp == "ols12 det(mult)") theta = (
		/* t_0    0.01 */  -3.9559128,   2.1481761,  -5.6798976,   6.2843492\   
		/*        0.05 */  -3.4087024,   2.1353073,  -2.8526172,   10.533851\   
		/*        0.10 */  -3.1262605,   2.0826855,  -1.8548174,   11.965616\   
		/* t_S/2  0.01 */  -3.9600056,   2.3310207,  -8.2349019,   17.322549\   
		/*        0.05 */  -3.4120999,   2.2749461,  -4.6162299,   17.309801\   
		/*        0.10 */  -3.1281091,   2.1448211,  -2.4136137,   12.862881\   
		/* F_K    0.01 */   12.214080,  -14.953704,   37.803806,  -35.469902\   
		/*        0.05 */   9.7501738,  -13.736069,   21.923671,  -48.117166\   
		/*        0.10 */   8.5753322,  -12.681100,   14.044172,  -38.983239\   
		/* F_SEAS 0.01 */   7.9955643,   5.5364164,   9.6733272,   66.803089\   
		/*        0.05 */   7.1481837,   2.2480373,   3.1483006,   14.525138\   
		/*        0.10 */   6.7192513,  0.80714794,  0.59725141,  -1.3134692\   
		/* F_ALL  0.01 */   7.8673893,   6.5760940,   10.458514,   64.373974\   
		/*        0.05 */   7.0632075,   3.0625000,   4.7586321,   7.0607786\   
		/*        0.10 */   6.6546712,   1.5835817,  0.27545972,   3.4407029    
	)
	/* CASE 1 */ if (tmp == "gls12 det(const)") theta = (
		/* t_0    0.01 */  -2.6064538,  -10.175621,   107.06159,  -402.02016\   
		/*        0.05 */  -1.9928579,  -12.360973,   127.92231,  -477.93486\   
		/*        0.10 */  -1.6750552,  -13.950681,   143.98809,  -538.39363\   
		/* t_S/2  0.01 */  -2.5653690,   1.0884361,  -2.5664519,   12.626518\   
		/*        0.05 */  -1.9399909,  0.88741346, -0.83355204,   5.3743644\   
		/*        0.10 */  -1.6159129,  0.77305443, -0.56572100,   4.4067282\   
		/* F_K    0.01 */   4.7332408,  -2.6292837,  -1.6384555,   15.680200\   
		/*        0.05 */   3.1100621,  -2.0903709,  -2.0046105,   16.114836\   
		/*        0.10 */   2.4080405,  -1.7772591, -0.77557285,   7.3328665\   
		/* F_SEAS 0.01 */   2.3426785,  0.42808389,  -4.0327865,   23.752202\   
		/*        0.05 */   1.8781229, -0.19242857, -0.86833969,   5.3513218\   
		/*        0.10 */   1.6547904, -0.28778292,  -1.6601140,   8.7904756\   
		/* F_ALL  0.01 */   2.2801347,   3.5575784,  -20.034165,   55.871893\   
		/*        0.05 */   1.8413546,   2.8052868,  -18.320792,   46.396554\   
		/*        0.10 */   1.6302099,   2.5580081,  -18.336290,   46.626318    
	)
	/* CASE 2 */ if (tmp == "gls12 det(trend)") theta = (
		/* t_0    0.01 */  -3.4304315,  -7.4282461,   91.886777,  -356.69718\   
		/*        0.05 */  -2.8781125,  -8.4054363,   104.26608,  -402.18423\   
		/*        0.10 */  -2.5934679,  -9.1778755,   114.42316,  -445.30801\   
		/* t_S/2  0.01 */  -2.5646271,  0.77224587,  0.39600462,  -2.0871711\   
		/*        0.05 */  -1.9396988,  0.65851311,   1.0901346,  -4.1652911\   
		/*        0.10 */  -1.6155446,  0.57221896,   1.1731581,  -4.6940397\   
		/* F_K    0.01 */   4.7325103,  -2.8332203,   2.1462070,   6.5769371\   
		/*        0.05 */   3.1102995,  -2.2223455,  0.39148271,   7.1689207\   
		/*        0.10 */   2.4073859,  -1.7890958,  -1.0971038,   13.233966\   
		/* F_SEAS 0.01 */   2.3441970,  0.24263481, -0.16043800,   7.1934005\   
		/*        0.05 */   1.8770934, -0.16308213,  -1.4027531,   11.669285\   
		/*        0.10 */   1.6550325, -0.34899254, -0.53525570,   5.8125377\   
		/* F_ALL  0.01 */   2.5458633,   5.4555905,  -48.924177,   195.25074\   
		/*        0.05 */   2.0850750,   4.6714437,  -49.067748,   194.37183\   
		/*        0.10 */   1.8627893,   4.3530181,  -48.705273,   192.25109    
	)
	/* CASE 3 */ if (tmp == "gls12 det(seas)") theta = (
		/* t_0    0.01 */  -2.6032254,  -10.565431,   89.999290,  -324.73774\   
		/*        0.05 */  -1.9900710,  -12.720026,   114.68568,  -416.39375\   
		/*        0.10 */  -1.6723983,  -14.264588,   132.26462,  -483.87790\   
		/* t_S/2  0.01 */  -2.6033922,  -10.571855,   90.112070,  -326.24338\   
		/*        0.05 */  -1.9882113,  -12.875656,   117.75139,  -433.64868\   
		/*        0.10 */  -1.6700677,  -14.477611,   136.40506,  -506.43205\   
		/* F_K    0.01 */   4.7221492,   22.174520,  -129.90660,   361.83506\   
		/*        0.05 */   3.0965730,   20.490423,  -109.26157,   251.84211\   
		/*        0.10 */   2.3935497,   19.062754,  -90.431116,   167.21374\   
		/* F_SEAS 0.01 */   2.3262363,   19.097837,  -44.878795,   56.146641\   
		/*        0.05 */   1.8606476,   16.958423,  -38.289220,   16.810367\   
		/*        0.10 */   1.6390210,   15.815438,  -33.022949,  -9.2209065\   
		/* F_ALL  0.01 */   2.2664664,   20.290791,  -51.106630,   72.442815\   
		/*        0.05 */   1.8262331,   18.238574,  -47.567433,   47.302576\   
		/*        0.10 */   1.6157800,   17.136718,  -43.799072,   28.997292    
	)
	/* CASE 4 */ if (tmp == "gls12 det(strend)") theta = (
		/* t_0    0.01 */  -3.4262229,  -7.9358366,   72.645255,  -264.86438\   
		/*        0.05 */  -2.8737864,  -8.9013045,   89.260973,  -338.73215\   
		/*        0.10 */  -2.5904240,  -9.5647527,   99.205847,  -380.26877\   
		/* t_S/2  0.01 */  -2.6033096,  -10.851371,   91.946728,  -337.23634\   
		/*        0.05 */  -1.9890582,  -12.985181,   116.53786,  -426.26983\   
		/*        0.10 */  -1.6704546,  -14.597529,   135.88504,  -502.82602\   
		/* F_K    0.01 */   4.7219074,   22.132648,  -126.63915,   356.32214\   
		/*        0.05 */   3.0964846,   20.463844,  -106.47170,   239.52831\   
		/*        0.10 */   2.3940686,   19.002652,  -86.878466,   147.08193\   
		/* F_SEAS 0.01 */   2.3260625,   19.106758,  -42.357516,   44.817087\   
		/*        0.05 */   1.8601860,   17.015625,  -37.070350,   11.463137\   
		/*        0.10 */   1.6382999,   15.895696,  -32.409380,  -12.082042\   
		/* F_ALL  0.01 */   2.5298922,   22.052039,  -73.866840,   189.47025\   
		/*        0.05 */   2.0693866,   19.910300,  -71.399675,   168.54082\   
		/*        0.10 */   1.8476919,   18.791946,  -68.234817,   151.99540    
	)
	/* CASE 5 */ if (tmp == "gls12 det(mult)") theta = (
		/* t_0    0.01 */  -3.4255075,  -10.916232,   78.735441,  -301.87954\   
		/*        0.05 */  -2.8735277,  -11.421414,   93.029268,  -359.18202\   
		/*        0.10 */  -2.5897606,  -11.918101,   103.72760,  -405.92011\   
		/* t_S/2  0.01 */  -3.4306407,  -10.619913,   73.569174,  -276.34133\   
		/*        0.05 */  -2.8764177,  -11.285219,   90.881955,  -348.71474\   
		/*        0.10 */  -2.5927494,  -11.781383,   101.68756,  -396.43024\   
		/* F_K    0.01 */   8.6680306,   41.175588,  -206.52703,   763.13099\   
		/*        0.05 */   6.5991203,   39.046128,  -242.56219,   878.43064\   
		/*        0.10 */   5.6375963,   38.436568,  -261.67489,   952.63835\   
		/* F_SEAS 0.01 */   5.2745442,   46.189689,  -199.46349,   862.48093\   
		/*        0.05 */   4.6063816,   43.236638,  -217.73805,   876.83087\   
		/*        0.10 */   4.2719648,   41.957602,  -226.62313,   887.03499\   
		/* F_ALL  0.01 */   5.1863965,   47.703900,  -211.37367,   928.73124\   
		/*        0.05 */   4.5515347,   44.753360,  -231.17955,   948.01106\   
		/*        0.10 */   4.2336770,   43.432724,  -240.63387,   962.66299    
	)

	// return critical values                                                   
	return(colshape(theta[,1] + theta[,2] * (N/s)^-1 + theta[,3] * (N/s)^-2    +
	                theta[,4] * (N/s)^-3, 3))
}
end

version 15.1: lmbuild larimaauto.mlib, replace size(8)
