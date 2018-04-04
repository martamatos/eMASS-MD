(* Created with the Wolfram Language for Students - Personal Use Only : www.wolfram.com *)
(2*metabolite["g3p", "c"]*metabolite["s7p", "c"]*parameter["TALA2_total"]*
  rateconst["TALA21", True]*rateconst["TALA22", True]*
  rateconst["TALA23", True]*rateconst["TALA24", True]*
  rateconst["TALA25", True]*rateconst["TALA26", True]*
  rateconst["TALA2_Kic_pi_1_e4p", False]*rateconst["TALA2_Kic_pi_2_e4p", 
   False])/(metabolite["s7p", "c"]*rateconst["TALA21", True]*
   rateconst["TALA24", True]*(rateconst["TALA23", False]*
     rateconst["TALA25", False] + rateconst["TALA22", True]*
     (rateconst["TALA23", False] + rateconst["TALA25", True]))*
   rateconst["TALA26", True]*(rateconst["TALA2_Kic_pi_1_e4p", False] + 
    metabolite["pi", "c"]*rateconst["TALA2_Kic_pi_1_e4p", True])*
   rateconst["TALA2_Kic_pi_2_e4p", False] + metabolite["g3p", "c"]*
   rateconst["TALA23", True]*rateconst["TALA2_Kic_pi_1_e4p", False]*
   (metabolite["s7p", "c"]*rateconst["TALA21", True]*
     (rateconst["TALA24", True]*(rateconst["TALA25", False] + 
        rateconst["TALA25", True])*rateconst["TALA26", True] + 
      rateconst["TALA22", True]*(rateconst["TALA24", False]*
         rateconst["TALA25", True] + rateconst["TALA25", True]*
         rateconst["TALA26", True] + rateconst["TALA24", True]*
         (rateconst["TALA25", True] + rateconst["TALA26", True])))*
     rateconst["TALA2_Kic_pi_2_e4p", False] + rateconst["TALA22", True]*
     rateconst["TALA24", True]*rateconst["TALA25", True]*
     rateconst["TALA26", True]*(rateconst["TALA2_Kic_pi_2_e4p", False] + 
      metabolite["pi", "c"]*rateconst["TALA2_Kic_pi_2_e4p", True]) + 
    rateconst["TALA21", False]*rateconst["TALA22", True]*
     rateconst["TALA25", True]*(rateconst["TALA24", False] + 
      rateconst["TALA26", True])*(rateconst["TALA2_Kic_pi_2_e4p", False] + 
      metabolite["pi", "c"]*rateconst["TALA2_Kic_pi_2_e4p", True])))
