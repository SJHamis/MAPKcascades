function monotherapy_dose = generate_monotherapy_dose_data(BRAF_in, ATP_in, IC_in, tend)

%Add a path to the directory that contains the system of reactions (model 
% structure) and model parameters.
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);

M = eye(46);
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
y0(42)=phosph3_tot;%phosph3
y0(15)=DBF_tot; %dbf
y0(30)=TMT_tot; %tmt
y0(44)=SCH_tot; %sch772984

IC_in = IC_in/100;

options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-3,'AbsTol',1e-3);%,'Vectorized','on');

monotherapy_points = zeros(3);


% Generate monotherapy doses for DBF, TMT, and SCH772984

% TMT monotherapy 
TMT_in = TMT_tot;
stage = 1;  % used to speed up code

while monotherapy_points(1,1) == 0
    y0(30) = TMT_in;
    [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_in, ATP_in, DBF_tot, TMT_in, SCH_tot), tspan, y0, options);
    ysubact = y(:,41)/1.2;
    t_index = find(t==tend*3600);
    ysubact = ysubact(t_index);

    if ysubact > IC_in && stage==1
        previous_ysubact = ysubact;
        previous_TMT = TMT_in;
        TMT_in = TMT_in+1;
    elseif ysubact <= IC_in && stage==1
        TMT_in=previous_TMT;
        stage=2;

    elseif ysubact > IC_in && stage==2
        previous_ysubact = ysubact;
        previous_TMT = TMT_in;
        TMT_in = TMT_in+0.1;
    elseif ysubact <= IC_in && stage==2
        TMT_in=previous_TMT;
        stage=3;

    elseif ysubact > IC_in && stage==3
        previous_ysubact = ysubact;
        previous_TMT = TMT_in;
        TMT_in = TMT_in+0.01;
    elseif ysubact <= IC_in && stage==3
        if abs(previous_ysubact-IC_in) <= abs(ysubact-IC_in)
            monotherapy_points(1,1) = previous_TMT;
            y0(30) = TMT_tot;
        else
            monotherapy_points(1,1) = TMT_in;
            y0(30) = TMT_tot;
        end   
    end
end


% DBF monotherapy 
DBF_in = DBF_tot;
stage = 1;  % used to speed up code

while monotherapy_points(2,2) == 0
    y0(15) = DBF_in;
    [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_in, ATP_in, DBF_in, TMT_tot, SCH_tot), tspan, y0, options);
    ysubact = y(:,41)/1.2;
    t_index = find(t==tend*3600);
    ysubact = ysubact(t_index);

    if ysubact > IC_in && stage==1
        previous_ysubact = ysubact;
        previous_DBF = DBF_in;
        DBF_in = DBF_in+1;
    elseif ysubact <= IC_in && stage==1
        DBF_in=previous_DBF;
        stage=2;

    elseif ysubact > IC_in && stage==2
        previous_ysubact = ysubact;
        previous_DBF = DBF_in;
        DBF_in = DBF_in+0.1;
    elseif ysubact <= IC_in && stage==2
        DBF_in=previous_DBF;
        stage=3;

    elseif ysubact > IC_in && stage==3
        previous_ysubact = ysubact;
        previous_DBF = DBF_in;
        DBF_in = DBF_in+0.01;
    elseif ysubact <= IC_in && stage==3
        if abs(previous_ysubact-IC_in) <= abs(ysubact-IC_in)
            monotherapy_points(2,2) = previous_DBF;
            y0(15) = DBF_tot;
        else
            monotherapy_points(2,2) = DBF_in;
            y0(15) = DBF_tot;
        end   
    end
end


% SCH772984 monotherapy 
SCH_in = SCH_tot;
stage = 1;  % used to speed up code

while monotherapy_points(3,3) == 0
    y0(44) = SCH_in;
    [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_in, ATP_in, DBF_tot, TMT_tot, SCH_in), tspan, y0, options);
    ysubact = y(:,41)/1.2;
    t_index = find(t==tend*3600);
    ysubact = ysubact(t_index);

    if ysubact > IC_in && stage==1
        previous_ysubact = ysubact;
        previous_SCH = SCH_in;
        SCH_in = SCH_in+1;
    elseif ysubact <= IC_in && stage==1
        SCH_in=previous_SCH;
        stage=2;

    elseif ysubact > IC_in && stage==2
        previous_ysubact = ysubact;
        previous_SCH = SCH_in;
        SCH_in = SCH_in+0.1;
    elseif ysubact <= IC_in && stage==2
        SCH_in=previous_SCH;
        stage=3;

    elseif ysubact > IC_in && stage==3
        previous_ysubact = ysubact;
        previous_SCH = SCH_in;
        SCH_in = SCH_in+0.01;
    elseif ysubact <= IC_in && stage==3
        if abs(previous_ysubact-IC_in) <= abs(ysubact-IC_in)
            monotherapy_points(3,3) = previous_SCH;
            y0(44) = SCH_tot;
        else
            monotherapy_points(3,3) = SCH_in;
            y0(44) = SCH_tot;
        end   
    end
end

monotherapy_dose = monotherapy_points;

end