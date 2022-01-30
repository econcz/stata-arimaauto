version 13.0

// CERTIFICATION FILE for ARIMAAuto()
discard
clear all
mata: AA = ARIMAAuto()

// quarterly data
sysuse gnp96.dta, clear

/* put variables */
mata: AA.put("varlist","gnp96")
mata: AA.put("trace", 2)
mata: AA.put("MS", ("",""))
/* intermediary results */
mata: AA.get("L")
mata: AA.get("T")
mata: AA.get("MS")
/* estimate */
mata: AA.start()
/* final results */
mata: AA.get("L")
mata: AA.get("T")
mata: AA.get("MS")

// monthly data
webuse air2.dta, clear

/* put variables */
mata: AA.put("varlist","air")
mata: AA.put("ifin","in 1/100")
mata: AA.put("MS", ("","0 1 0 12"))
/* intermediary results */
mata: AA.get("L")
mata: AA.get("T")
mata: AA.get("MS")
/* estimate */
mata: AA.start()
/* final results */
mata: AA.get("L")
mata: AA.get("T")
mata: AA.get("MS")

// CLEAR MEMORY
mata: mata drop AA
