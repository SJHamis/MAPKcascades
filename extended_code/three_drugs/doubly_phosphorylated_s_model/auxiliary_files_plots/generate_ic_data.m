function IC_value = generate_ic_data(BRAF_in, ATP_in, IC_in, tend)

%Add a path to the directory that contains the system of reactions (model 
% structure) and model parameters.
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);

M = eye(51);
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
SUB_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
phosph3_tot=0.12;
DBF_tot=0;
TMT_tot=0;
SCH_tot=0;

y0(1)=BRAF_tot;%braf
y0(2)=ATP_tot;%atp
y0(4)=MEK_tot;%mek
y0(20)=ERK_tot; %erk
y0(38)=SUB_tot;%substrate sub
y0(12)=phosph1_tot;%phosph1
y0(27)=phosph2_tot;%phosph2
y0(45)=phosph3_tot;%phosph3
y0(15)=DBF_tot; %dbf
y0(30)=TMT_tot; %tmt
y0(48)=SCH_tot; %sch772984

options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-3,'AbsTol',1e-3);%,'Vectorized','on');

IC_in = IC_in/100;
[t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_in, ATP_in, DBF_tot, DBF_tot, DBF_tot), tspan, y0, options);
ypsub = y(:,40)/1.2;
yppsub = y(:,44)/1.2;
ysubact_nodrugs = ypsub+yppsub;
t_index = find(t==tend*3600);
ysubact_nodrugs = ysubact_nodrugs(t_index);
IC_value = IC_in*ysubact_nodrugs;

end