(* Created with the Wolfram Language for Students - Personal Use Only : www.wolfram.com *)
(metabolite["2pg", "c"]*rateconst["ENO1", True]*
  (rateconst["ENO2", True] + rateconst["ENO3", False] + 
   rateconst["ENO3", True]))/
 (rateconst["ENO1", False]*(rateconst["ENO2", True] + 
    rateconst["ENO3", False]) + rateconst["ENO2", True]*
   rateconst["ENO3", True] + metabolite["2pg", "c"]*rateconst["ENO1", True]*
   (rateconst["ENO2", True] + rateconst["ENO3", False] + 
    rateconst["ENO3", True]))
