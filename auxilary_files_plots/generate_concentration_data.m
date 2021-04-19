function pp_fraction_ss = generate_concentration_data(BRAF_in, ATP_in, DBF_in, TMT_in, tend)

%Add a path to the directory that contains the system of reactions (model 
% structure) and model parameters.
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);

M = eye(36);
super_compound_list_index=get_conslaw_position();

%Substitute in cons. laws 
for i=1:size(M,1)
    if(ismember(i,super_compound_list_index))
        M(i,i)=0; 
    end
end

tspan = [0 tend*3600];
y0=zeros(size(M,1),1);

BRAF_tot=BRAF_in;
ATP_tot=ATP_in;
MEK_tot=1.2;
ERK_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
DBF_tot=DBF_in;
TMT_tot=TMT_in;


y0(1)=BRAF_tot;%braf
y0(2)=ATP_tot;%atp
y0(4)=MEK_tot;%mek
y0(20)=ERK_tot; %erk
y0(12)=phosph1_tot;%phosph1
y0(27)=phosph2_tot;%phos2
y0(15)=DBF_tot; %dbf
y0(30)=TMT_tot; %tmt

options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-3,'AbsTol',1e-3);%,'Vectorized','on');

[~,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_in, ATP_in, DBF_in, TMT_in), tspan, y0, options);

%Doubly phosphorylated fraction at steady state: ;
pp_fraction_ss=[y(end,26)/ERK_tot];

end