function [system_of_ODEs_Math_LHS_out,system_of_ODEs_Math_RHS_out,system_of_ODEs_text_out] = generate_odes_from_reactions(list_of_compounds_in,srmf_math_in,srmf_text_in)

% Generate an output system of ODEs from the input system of reactions 
list_of_compounds=list_of_compounds_in;
srmf_math=srmf_math_in;
srmf_text=srmf_text_in;
system_of_odes_LHS=[];
system_of_odes_RHS=[];
system_of_odes_text=[];

%srmf_math/text(z,r,c) returns the r:th row and the c:th column entry of
%srmf_math/text
z=1;
norows=size(srmf_math,2);
%nocols=size(srmf_math,3);


for compound_i = 1:length(list_of_compounds)
    %for every row, see if the i:th compund occurrs in cols 1,4,6 i.e
    %reactansts(co11). complex(col4) product (col6)
    temp_string = "";
    temp_string_text = "";
    compound_i_text = regexprep(list_of_compounds(compound_i), '\\[', '[');
    compound_i_text = regexprep(compound_i_text, '\\]', ']');
    
    for r = 1:norows
        reactants=srmf_math(z,r,1);
        a=srmf_math(z,r,2);
        d=srmf_math(z,r,3);
        complex=srmf_math(z,r,4);
        k=srmf_math(z,r,5);
        product=srmf_math(z,r,6);
        
        reactants_text=srmf_text(z,r,1);
        a_text=srmf_text(z,r,2);
        d_text=srmf_text(z,r,3);
        complex_text=srmf_text(z,r,4);
        k_text=srmf_text(z,r,5);
        product_text=srmf_text(z,r,6);
       
        %if compound i in list_of_compounds exists in reactants 
        if(regexp(reactants,'.*y\('+string(compound_i)+'\).*'))
             % add the term: 
             %(product of reactants)*(-col2) + col(4)*col(3)
             %where col 2 is a-value, col3 is d-value and col4 is complex
             %(product of reactants)*(-srmf_math(z,r,2)) + srmf_math(z,r,3)*srmf_math(z,r,4)
             temp_reactant_product=regexprep(reactants,"+","*");
             %temp_string = temp_string + "+" + "y(" + string(compound_i) + ")" + "*" + "(-" + a + ")+" + d + "*" + complex;
             temp_string = temp_string + "+" + temp_reactant_product + "*" + "(-" + a + ")+" + d + "*" + complex;
        end

        %if compound i in list_of_compounds exists in complex
        if(regexp(complex,'.*y\('+string(compound_i)+'\).*'))
             % add the term:
             %y(i)*(-col3-col5) + (col2)*prod(col1) 
             %where col 3 is d-value, col5 is k-value, col2 is a-value and col1 is reactants
             %prod(col1) is the products of all complexes in the reactants
             temp_string = temp_string + "+" + "y(" + string(compound_i) +")" + "*" + "(" + "-" + d + "-" + k + ")"...
                 + "+" + a + "*" + "(" + regexprep(reactants,"+","*") + ")"; 
             %where regexprep(srmf_math(z,r,1),"+","*") yields the product of reactanst
        end

        %if compound i in list_of_compounds exists in product
        if(regexp(product,'.*y\('+string(compound_i)+'\).*'))
             % add the term:
             % col4*col5
             % where col4 is the complex and col5 is the k-value
             temp_string = temp_string + "+" + complex + "*" + k;
        end
        
        %%%%%%%%%%%%%%%%%
        %%%text version
        %%%%%%%%%%%%%%%%%
        
        %if compound i in list_of_compounds exists in reactants 
        if(contains(reactants_text,compound_i_text))
             temp_reactant_product_text=regexprep(reactants_text,"+","*");
             temp_string_text = temp_string_text + "+" + temp_reactant_product_text + "*" + "(-" + a + ")+" + d + "*" + complex_text;
       end

        %if compound i in list_of_compounds exists in complex
        if(contains(complex_text,compound_i_text))
             temp_string_text = temp_string_text + "+" + compound_i_text + "*" + "(" + "-" + d + "-" + k + ")"...
                 + "+" + a + "*" + "(" + regexprep(reactants_text,"+","*") + ")"; 
        end

        %if compound i in list_of_compounds exists in product
         if(contains(product_text,compound_i_text))
             temp_string_text = temp_string_text + "+" + complex_text + "*" + k;
        end
    end
    
    %add the temporary stings into systems
    temp_string_LHS = "dydt(" + string(compound_i) + ")";
    system_of_odes_LHS=[system_of_odes_LHS; temp_string_LHS];
    system_of_odes_RHS=[system_of_odes_RHS; temp_string];
    temp_string_text = "d" + compound_i_text +"/dt =" + temp_string_text;
    system_of_odes_text=[system_of_odes_text; temp_string_text];
end

%remove -k0 dummy entries from system of ODEs 
system_of_odes_RHS=regexprep(system_of_odes_RHS,"-k0","");
system_of_odes_text=regexprep(system_of_odes_text,"-k0","");
system_of_odes_text=regexprep(system_of_odes_text,"=\+","= ");

system_of_ODEs_Math_LHS_out=system_of_odes_LHS;
system_of_ODEs_Math_RHS_out=system_of_odes_RHS;
system_of_ODEs_text_out=system_of_odes_text;


end

