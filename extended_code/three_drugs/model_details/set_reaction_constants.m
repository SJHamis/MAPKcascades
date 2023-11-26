function [reaction_constants] =  set_reaction_constants()%

total_concentrations = ["BRAF_tot=BRAF_in;"
    "ATP_tot=ATP_in;"
    "MEK_tot=1.2;"
    "ERK_tot=1.2;"
    "SUB_tot=1.2;"   % SUB = dummy substrate
    "phosph1_tot=0.0003;"
    "phosph2_tot=0.12;"
    "phosph3_tot=0.12;"  
    "DBF_tot=DBF_in;"
    "TMT_tot=TMT_in;"
    "SCH_tot=SCH_in;"];

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
    "a10=a_HF;"
    "a11=a_HF;"
    "a12=a_HF;"];   

d_backward_constants = ["d1=6.23; %0.8268; % ATP binding to BRAF (or BRAF complex)."
"d2=0.02385; % MEK binding to BRAF (or BRAF complex)."
"d5=0.3468; % ATP binding to MEK (or MEK complex)."
"d6=0.01184; % ERK binding to MEK (or MEK complex)."
"d4=0.00005936; % DBF binding to BRAF (or BRAF complex)."
"d8=0.0012296; % TMT binding to MEK (or MEK complex)."
"d3=0.0159; % BRAF phosphasate."
"d7=0.0159; % MEK phosphasate."
"d9=14.64; % ATP binding to ERK (or ERK complex)."
"d10=0.7738; % SUB binding to ERK (or ERK complex)."
"d11=0.0159; % ERK phosphatase."
"d12=0.0002967; % SCH binding to ERK (or ERK complex)."];  

k_cat_constants = ["k12=0.66;"
"k56=0.0242;"
"k3=0.0159;"
"k7=0.0159;"
"k910=0.2;"
"k11=0.0159;"];   
   
reaction_constants=[total_concentrations; a_forward_constants;...
    d_backward_constants; k_cat_constants];

end