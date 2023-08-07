function [reaction_constants] =  set_reaction_constants()%

total_concentrations = ["KKK_tot=KKK_in;"
    "E1_tot=E1_in;"
    "E2_tot=E2_in;"
    "KK_tot=1.2;"
    "K_tot=1.2;"
    "phosph1_tot=0.0003;"
    "phosph2_tot=0.12;"];

a_forward_constants = ["a_HF=0.106; %All a-values are the same.";
    "a1=a_HF;"
    "a2=a_HF;"
    "a3=a_HF;"
    "a4=a_HF;"
    "a5=a_HF;"
    "a6=a_HF;"
    "a7=a_HF;"
    "a8=a_HF;"
    "a9=a_HF;"
    "a10=a_HF;"];

d_backward_constants = ["d1=6.23; %0.8268; % E1 and E2 binding to KKK (or KKK complex)."
"d2=d1;"
"d3=0.02385; % KK binding to KKK (or KKK complex)."
"d5=d3;"
"d7=0.01184; % K binding to KK (or KK complex)."
"d9=d7;"
"d4=0.0159; % KKK phosphasate."
"d6=d4;"
"d8=0.0159; % KK phosphasate."
"d10=d8;"]; 

k_cat_constants = ["k1=0.66;"
"k2=k1;"
"k3=k1;"
"k5=k1;"
"k7=0.0242;"
"k9=k7;"
"k4=0.0159;"
"k6=k4;"
"k8=0.0159;"
"k10=k8;"]; 
   
reaction_constants=[total_concentrations; a_forward_constants;...
    d_backward_constants; k_cat_constants];

end