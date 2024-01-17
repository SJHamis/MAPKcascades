close all
clear all

%Add a path to the directory that contains the model details
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

% Set inital conditions
BRAF_tot=0.003;
MEK_tot=1.2;
ERK_tot=1.2;
SUB_tot=1.2;
ATP_tot=1000;
phosph1_tot=0.0003;
phosph2_tot=0.12;
phosph3_tot=0.12;
DBF_tot=0;
TMT_tot=0;

y0=zeros(size(M,1),1);
y0(1)=BRAF_tot;
y0(4)=MEK_tot;
y0(20)=ERK_tot;
y0(38)=SUB_tot;
y0(2)=ATP_tot;
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;
y0(42)=phosph3_tot;
y0(15)=DBF_tot;
y0(30)=TMT_tot;


tspan = [0 3*3600];   % 3h (can be changed later)
SCH_tot=[0 1 5 10 16.2 25 50 100]/10;

options = odeset('Mass',M, 'MassSingular','yes', 'RelTol',1e-3, 'AbsTol',1e-3);
set(gcf,'position',[100,100,1200,700])
newcolors = {'#0000FF','#008000',' #FFC300 ', '#FF5733 ','#C70039','#900C3F','#581845','#000000','#DAF7A6'};

%%%%%%%%%%%%%%%%%%%%%%% Generate activated substrate SUB (unbound pSUB) vs t FIGURE for different SCH772984 values

colororder(newcolors)

for j = 1:length(SCH_tot)
    SCH_in=SCH_tot(j);
    y0(44)=SCH_in;
    [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot, ATP_tot, DBF_tot, TMT_tot, SCH_in), tspan, y0, options);
    ysubact = y(:,41)/1.2;
    hold on
    plot(t,ysubact,'LineWidth',2.0)
    clear t y
end

xticks(3600*[0 1 2 3])
xlim([0 tspan(end)])
xticklabels({'0','1','2','3'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',12)
ylim([0 1])
grid on
grid minor

lgd = legend('0\muM', '0.1\muM', '0.5\muM', '1\muM', '1.62\muM', '2.5\muM', '5\muM', '10\muM', ...
'Orientation', 'horizontal', 'Location','southoutside');
title(lgd, "SCH doses");