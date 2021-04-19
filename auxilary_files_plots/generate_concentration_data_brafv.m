function reqBRAFlevel = generate_concentration_data_brafv(DBF_in, TMT_in, tend)

%Get current directory location
[parentdir, ~,~]=fileparts(pwd);
%Directory of model files 
dae_location = strcat(parentdir,'/auxiliary_files_model_setup');
addpath(dae_location);

M = eye(36);
super_compound_list_index=get_conslaw_position;


%Substitute in cons. laws 
for i=1:size(M,1)
    if(ismember(i,super_compound_list_index))
        M(i,i)=0; 
    end
end

tspan = [0 tend*3600];
y0=zeros(size(M,1),1);


ATP_tot=1000;
MEK_tot=1.2;
ERK_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
DBF_tot=DBF_in;
TMT_tot=TMT_in;



y0(2)=ATP_tot;%atp
y0(4)=MEK_tot;%mek
y0(20)=ERK_tot; %erk
y0(12)=phosph1_tot;%phosph1
y0(27)=phosph2_tot;%phos2
y0(15)=DBF_tot; %dbf
y0(30)=TMT_tot; %tmt

options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-3,'AbsTol',1e-3);%,'Vectorized','on');
min=0.1*0.001;
max=30*0.001;
points=3000;
BRAF_values = linspace(min,max,points);
BRAF_index = 1;
goal_value_reached = 0;

while(goal_value_reached==0 && BRAF_index<points)
    BRAF_in = BRAF_values(BRAF_index);
    y0(1)=BRAF_in; 
    [~,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_in, ATP_tot, DBF_in, TMT_in), tspan, y0, options);
    pp_fraction_ss=[y(end,26)/1.2];
    if(pp_fraction_ss>0.5)
        reqBRAFlevel = BRAF_in;
        goal_value_reached=1;
    else
        BRAF_index = BRAF_index+1;
    end
end

if(BRAF_index>=points)
    reqBRAFlevel=NaN;%BRAF_values(end);
end


end