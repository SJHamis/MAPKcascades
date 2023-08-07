clear all
close all

%Add a path to the directory that contains the system of reactions (model 
% structure) and model parameters.
%Get current directory location
[parentdir, ~,~]=fileparts(pwd);
%Directory of model files 
model_location = strcat(parentdir,'/model_details');
addpath(model_location);


% Step 1: Enter the system of reactions (SoRs) manually in the 
% MATLAB-function-file  set_system_of_reactions.m.

% Step 2: Convert the SoRs into a system of ordinary differential equations
% (SoODEs) by running the (this) MATLAB file 
% (main_setup_convert_SoRs_to_SoODEs.m).

% Read in a 1xN string (where N = the number of dependent variables
% in the system of reactions) from the set_system_of_reaction.m. 
system_of_reactions = set_system_of_reactions();

% Here, system_of_reactions(i) corresponds to the i:th reaction. 
% The reactions are written in form [A]+[B] aj dj [A_B] kj [A]+[C],
% where k0 and [empty] denotes zero valued/non existing components. 
% [A] denotes the concentration of A. [A_B] denotes the concentration of
% the complex of A and B.

% Using the MATLAB-function-file sort_complex_names, 
% we re-write the reaction equations so that every complex is in an 
% alphabetically sorted form (in order to avoid with duplicates). 
% For example: [B_A] becomes [A_B], so that [B_A] and [A_B] are
% recognised as the same complex.
for row=1:size(system_of_reactions,2)
    system_of_reactions(row)=sort_complex_names(system_of_reactions(row));
end

% We now create a list of all compounds (list_of_compounds) in the SoRs.
% We start with an empty string:
list_of_compounds="";

% Now search the SoRs for all compounds (occurences in square brackets []), 
% and add these compounds to the list.
for i = 1:size(system_of_reactions,2)
    temp_string = extractBetween(system_of_reactions(i),"[","]");
    list_of_compounds = [list_of_compounds;"\["+temp_string+"\]"];
end 

% Now remove all non-unique list entries using the function unique.
list_of_compounds=unique(list_of_compounds,'stable');
% The 'stable' option ensures that the compounds are listed in order of 
% apperance (not default alphabetically).

% Remove list entries that are named "empty" or are empty
list_of_compounds=erase(list_of_compounds,'\[empty\]');
list_of_compounds(cellfun('isempty',list_of_compounds))=[];

% Okay, so now we have the SoRs and a list of the N unique compounds.
% We re-write the i:th compound in the list (list_of_compounds) as y(i), so
% that we can turn the SoRs into a form suitable for MATLAB's ode-solvers.
% This gives us a modified reaction system in "math form", called 
% system_of_reactions_math_form.

% Initial form of system_of_reactions_math_form: 
system_of_reactions_math_form=system_of_reactions;

% Now update system_of_reactions_math_form... 
% For all compounds i in list_of_compounds
for i = 1:length(list_of_compounds)
    temp_string=list_of_compounds(i);
    i_string=string(i);
    % For all rows j in the SoRs
    for j = 1:size(system_of_reactions,2)
            %Change "[compound i in list_of_compounds]" to "y(i)" every
            %time "[compound i in list_of_compounds]" appears in the SoRs: 
            system_of_reactions_math_form(j) = ...
                regexprep(system_of_reactions_math_form(j), ...
                temp_string, "y("+i_string+")");
    end
end

% We split up system_of_reactions_math_form into matrix form for easy 
% sub-string manipulations:
srmf_math = split(system_of_reactions_math_form);
% Naming from (s)ystem of (r)eactions in (m)atrix (f)ormat (math) version.

% We also create a "text" version called
% (s)ystem of (r)eactions in (m)atrix (f)ormat (text) version, this is use-
% ful for readability, sanity checks and creating tex-formatted math 
% content.
srmf_text = split(system_of_reactions);

% We can now generate a SoODEs from the SoRs. 
% We get the LHS and RHS of the mathematical representation, as well as a
% text representation using the MATLAB-file-function
% generate_odes_from_reactions.
[system_of_ODEs_Math_LHS, system_of_ODEs_Math_RHS, system_of_ODEs_text]...
    = generate_odes_from_reactions(list_of_compounds,srmf_math,srmf_text);
% system_of_ODEs_Math_LHS looks like [dydt(1); dydt(2); ... dydt(N)]
% system_of_ODEs_Math_RHS may look like [+y(1)*y(2)*(-a1);...]
% Then system_of_ODEs_Math_text becomes [dydt(1)=+y(1)*y(2)*(-a1);...]

% We now incorporate the conservations laws associated with the SoRs.
% In out specific case, note that the total amount of 
% BRAF, MEK, ERK, O, ATP, phosph1, phosph2, phosph3, DBF, TMT, SCH (in 
% bound and unbound form) is conserved.
% For example, [BRAF]+[BRAF_A]+[BRAF_B] = constant = BRAF_tot,
% (or  d[BRAF]/dt + d[BRAF_A]/dt + d[BRAF_B]/dt = 0).
% If we re-write this as LHS = 0, we get
% [BRAF]+[BRAF_A]+[BRAF_B]-BRAF_tot = 0.

% We create a list of "super-compounds", containing every conservation-law
% related total-concentration-constant.
list_of_super_compounds=["BRAF";"MEK";"ERK";"O";"ATP";"DBF";...
    "TMT";"SCH";"phosph1";"phosph2";"phosph3"];

% Now find all sub-compounds to each super-compound by scanning through the
% list_of_compounds for occurenses of the super-compound using the MATLAB-
% function-file get_conservation_laws:
[conservation_laws_math_LHS, conservation_laws_math_RHS, ...
    conservation_laws_text, super_compound_list_index] = ...
    get_conservation_laws(list_of_super_compounds, list_of_compounds);

%super_compound_list_index=[1 2 4];

% We create an Nx1 string, which contains the RHS of the SoODEs. 
yp_tofile=system_of_ODEs_Math_RHS;

% We substitute in the conservation laws into the differential algebraic
% equation: 
j=1;
for i=super_compound_list_index
    yp_tofile(i)=conservation_laws_math_RHS(j);
    j=j+1;  
end


% Do some text formatting to be suitable for MATLAB-ode-solver (eg ode15s):
yp_tofile(1)="yp=["+yp_tofile(1);
yp_tofile(end)=yp_tofile(end)+"];";
% So we just added some formatting before and after
% system_of_ODEs_Math_RHS.

% We create a file called math_yp.m containing the yp_tofile string. 
file_name_math_yp = "mapk_cascade_DAE";
file_math_yp = fopen(file_name_math_yp+".m",'w');
% We add some formatting to this file: 
file_header_ode_math="function yp = " + file_name_math_yp + "(y, BRAF_in, ATP_in, DBF_in, TMT_in, SCH_in)";
file_TOT_constants_ode_math = set_reaction_constants();
file_end_ode_math="end";

fprintf(file_math_yp, '%s\n\n', file_header_ode_math);

for i = 1:size(file_TOT_constants_ode_math,1)
     fprintf(file_math_yp, '%s\n', file_TOT_constants_ode_math(i));
end
fprintf(file_math_yp, '\n');

for i = 1:size(yp_tofile,1)
     fprintf(file_math_yp, '%s\n', yp_tofile(i));
end
fprintf(file_math_yp, '\n%s', file_end_ode_math);
fclose(file_math_yp);


% Make a function that returns the positions in the result vector (y(t)) 
% where the conservation laws can be substituted in. 
file_name_conslaw_pos = "get_conslaw_position";
file_conslaw_pos = fopen(file_name_conslaw_pos +".m",'w');
% set strings
file_header_conslaw_pos="function clp = " + file_name_conslaw_pos + "()";
file_content_conslaw_pos = "clp = ["+num2str(super_compound_list_index)+"];";
file_end_conslaw_pos="end";
% write to file
fprintf(file_conslaw_pos, '%s\n\n', file_header_conslaw_pos);
fprintf(file_conslaw_pos, '%s\n\n', file_content_conslaw_pos);
fprintf(file_conslaw_pos, '\n%s', file_end_conslaw_pos);
fclose(file_conslaw_pos);





% Let's clear the workspace from unneccesary variables:
unnec_variables = {'row', 'i', 'temp_string', 'i_string', 'j', 'srmf_math', 'srmf_text', ...
'yp_tofile', 'super_compound_list_index', 'file_name_math_yp', 'file_math_yp'...
'file_header_ode_math', 'file_TOT_constants_ode_math', 'file_end_ode_math', 'ans' };
clear(unnec_variables{:})
clear unnec_variables

% We have now created a MATLAB-function file mapk_cascade_DAE.m with the
% cascade in ODE-form, with conservation laws. The DAE can be solved using eg. MATLAB's ODE solver 
% ode15s (as is done from out main files).