function [reactions_out] =  set_system_of_reactions()
% Set the chemical reactions.
% a: forward reaction
% b: backward reaction
% k: cat reaction
% If no cat reaction - write dummy term "k0 [empty]"

    reactions_KKK_cat = ["[KKK]+[E1] a1 d1 [KKK_E1] k1 [KKKact]+[E1];"...
        "[KKKact]+[E2] a2 d2 [KKK_E2] k2 [KKK]+[E2];"...
        "[KK]+[KKKact] a3 d3 [KK_KKKact] k3 [pKK]+[KKKact];"...
        "[pKK]+[KKKact] a5 d5 [pKK_KKKact] k5 [ppKK]+[KKKact];"...
        ]; 

    reactions_KKK_phosphotase = ["[pKK]+[phosph1] a4 d4 [pKK_phosph1] k4 [KK]+[phosph1];"...
        "[ppKK]+[phosph1] a6 d6 [ppKK_phosph1] k6 [pKK]+[phosph1]"...
        ]; 

    reactions_KK_cat = ["[ppKK]+[K] a7 d7 [ppKK_K] k7 [ppKK]+[pK];"...
        "[pK]+[ppKK] a9 d9 [pK_ppKK] k9 [ppK]+[ppKK];"...
        ]; 

    reactions_KK_phosphotase = ["[pK]+[phosph2] a8 d8 [pK_phosph2] k8 [K]+[phosph2];"...
        "[ppK]+[phosph2] a10 d10 [ppKK_phosph2] k10 [pK]+[phosph2]"...
        ]; 

    reactions_out=[reactions_KKK_cat, reactions_KKK_phosphotase,...
        reactions_KK_cat, reactions_KK_phosphotase];

end