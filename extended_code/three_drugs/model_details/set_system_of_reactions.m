function [reactions_out] =  set_system_of_reactions()
% Set the chemical reactions.
% a: forward reaction
% b: backward reaction
% k: cat reaction
% If no cat reaction - write dummy term "k0 [empty]"

    reactions_BRAF_cat = ["[BRAF]+[ATP] a1 d1 [BRAF_ATP] k0 [empty];"...
        "[BRAF]+[MEK] a2 d2 [BRAF_MEK] k0 [empty];"...
        "[BRAF]+[pMEK] a2 d2 [BRAF_pMEK] k0 [empty];"...
        "[BRAF_ATP]+[MEK] a2 d2 [BRAF_ATP_MEK] k12 [BRAF]+[ADP]+[pMEK];"...
        "[BRAF_ATP]+[pMEK] a2 d2 [BRAF_ATP_pMEK] k12 [BRAF]+[ADP]+[ppMEK];"...
        "[BRAF_MEK]+[ATP] a1 d1 [BRAF_MEK_ATP] k12 [BRAF]+[pMEK]+[ADP];"...
        "[BRAF_pMEK]+[ATP] a1 d1 [BRAF_pMEK_ATP] k12 [BRAF]+[ppMEK]+[ADP]"...
        ]; 

    reactions_BRAF_phosphotase = ["[pMEK]+[phosph1] a3 d3 [pMEK_phosph1] k3 [MEK]+[phosph1];"...
        "[ppMEK]+[phosph1] a3 d3 [ppMEK_phosph1] k3 [pMEK]+[phosph1]"...
        ]; 

    reactions_BRAF_inhibition = ["[BRAF]+[DBF] a4 d4 [BRAF_DBF] k0 [empty];"...
        "[BRAF_MEK]+[DBF] a4 d4 [BRAF_MEK_DBF] k0 [empty];"...
        "[BRAF_pMEK]+[DBF] a4 d4 [BRAF_pMEK_DBF] k0 [empty];"...
        "[BRAF_DBF]+[MEK] a2 d2 [BRAF_DBF_MEK] k0 [empty];"...
        "[BRAF_DBF]+[pMEK] a2 d2 [BRAF_DBF_pMEK] k0 [empty]"...
        ]; 

    reactions_MEK_cat = ["[ppMEK]+[ATP] a5 d5 [ppMEK_ATP] k0 [empty];"...
        "[ppMEK]+[ERK] a6 d6 [ppMEK_ERK] k0 [empty];"...
        "[ppMEK]+[pERK] a6 d6 [ppMEK_pERK] k0 [empty];"...
        "[ppMEK_ATP]+[ERK] a6 d6 [ppMEK_ATP_ERK] k56 [ppMEK]+[ADP]+[pERK];"...
        "[ppMEK_ATP]+[pERK] a6 d6 [ppMEK_ATP_pERK] k56 [ppMEK]+[ADP]+[ppERK];"...
        "[ppMEK_ERK]+[ATP] a5 d5 [ppMEK_ERK_ATP] k56 [ppMEK]+[pERK]+[ADP];"...
        "[ppMEK_pERK]+[ATP] a5 d5 [ppMEK_pERK_ATP] k56 [ppMEK]+[ppERK]+[ADP]"...
        ]; 

    reactions_MEK_phosphotase = ["[pERK]+[phosph2] a7 d7 [pERK_phosph2] k7 [ERK]+[phosph2];"...
        "[ppERK]+[phosph2] a7 d7 [ppERK_phosph2] k7 [pERK]+[phosph2]"...
        ]; 

    reactions_MEK_inhibition = ["[ppMEK]+[TMT] a8 d8 [ppMEK_TMT] k0 [empty];"...
        "[ppMEK_ATP]+[TMT] a8 d8 [ppMEK_ATP_TMT] k0 [empty];"...
        "[ppMEK_ERK]+[TMT] a8 d8 [ppMEK_ERK_TMT] k0 [empty];"...
        "[ppMEK_pERK]+[TMT] a8 d8 [ppMEK_pERK_TMT] k0 [empty];"...
        "[ppMEK_ERK_ATP]+[TMT] a8 d8 [ppMEK_ERK_ATP_TMT] k0 [empty];"...
        "[ppMEK_pERK_ATP]+[TMT] a8 d8 [ppMEK_pERK_ATP_TMT] k0 [empty];"...
        "[ppMEK_TMT]+[ATP] a5 d5 [ppMEK_TMT_ATP] k0 [empty];"...
        "[ppMEK_TMT]+[ERK] a6 d6 [ppMEK_TMT_ERK] k0 [empty];"...
        "[ppMEK_TMT]+[pERK] a6 d6 [ppMEK_TMT_pERK] k0 [empty];"... 
        "[ppMEK_TMT_ERK]+[ATP] a5 d5 [ppMEK_TMT_ERK_ATP] k0 [empty];"...
        "[ppMEK_TMT_pERK]+[ATP] a5 d5 [ppMEK_TMT_pERK_ATP] k0 [empty];"...
        "[ppMEK_TMT_ATP]+[ERK] a6 d6 [ppMEK_TMT_ATP_ERK] k0 [empty];"...
        "[ppMEK_TMT_ATP]+[pERK] a6 d6 [ppMEK_TMT_ATP_pERK] k0 [empty]"...
        ];

    reactions_ERK_cat = ["[ppERK]+[ATP] a9 d9 [ppERK_ATP] k0 [empty];"...
        "[ppERK]+[O] a10 d10 [ppERK_O] k0 [empty];"...
        "[ppERK]+[pO] a10 d10 [ppERK_pO] k0 [empty];"...
        "[ppERK_ATP]+[O] a10 d10 [ppERK_ATP_O] k910 [ppERK]+[ADP]+[pO];"...
        "[ppERK_ATP]+[pO] a10 d10 [ppERK_ATP_pO] k910 [ppERK]+[ADP]+[ppO];"...
        "[ppERK_O]+[ATP] a9 d9 [ppERK_ATP_O] k910 [ppERK]+[ADP]+[pO];"...
        "[ppERK_pO]+[ATP] a9 d9 [ppERK_ATP_pO] k910 [ppERK]+[ADP]+[ppO]"...
        ];

    reactions_ERK_phosphotase = ["[pO]+[phosph3] a11 d11 [pO_phosph3] k11 [O]+[phosph3];"...
        "[ppO]+[phosph3] a11 d11 [ppO_phosph3] k11 [pO]+[phosph3]"...
        ];

    reactions_ERK_inhibition = ["[ppERK]+[SCH] a12 d12 [ppERK_SCH] k0 [empty];"...
        "[ppERK_O]+[SCH] a12 d12 [ppERK_O_SCH] k0 [empty];"...
        "[ppERK_pO]+[SCH] a12 d12 [ppERK_pO_SCH] k0 [empty];"...
        "[ppERK_SCH]+[O] a10 d10 [ppERK_SCH_O] k0 [empty];"...
        "[ppERK_SCH]+[pO] a10 d10 [ppERK_SCH_pO] k0 [empty]"...
        ];

    reactions_out=[reactions_BRAF_cat, reactions_BRAF_phosphotase, reactions_BRAF_inhibition,...
        reactions_MEK_cat, reactions_MEK_phosphotase, reactions_MEK_inhibition,...
        reactions_ERK_cat, reactions_ERK_phosphotase, reactions_ERK_inhibition];

end