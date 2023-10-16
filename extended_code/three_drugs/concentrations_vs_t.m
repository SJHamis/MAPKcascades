close all
clear all

%Add a path to the directory that contains the model details
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

% Set inital conditions
MEK_tot=1.2;
ERK_tot=1.2;
O_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
phosph3_tot=0.12;
DBF_tot=0;
TMT_tot=0;
SCH_tot=2.5;

y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot;
y0(38)=O_tot;
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;
y0(45)=phosph3_tot;
y0(15)=DBF_tot;
y0(30)=TMT_tot;
y0(48)=SCH_tot;


tspan = [0 3*3600];   % 3h (can be changed later)

options = odeset('Mass',M, 'MassSingular','yes', 'RelTol',1e-3, 'AbsTol',1e-3);
set(gcf,'position',[100,100,1200,700])
newcolors = {'#FF00FF','#7E2F8E','#0000FF','#0072BD','#008000',' #FFC300 ', '#FF5733 ','#C70039','#900C3F','#581845','#000000','#DAF7A6'};

%%%%%%%%%%%%%%%%%%%%%%% Generate activated substrate O (pO + ppO) vs t FIGURE for different SCH772984 values
tiles = tiledlayout(2,2,'TileSpacing','compact'); 

BRAF_tot=[0.003,0.003,0.01,0.01];
ATP_tot=[1000, 5000, 1000, 5000];

for k = 1:4 
    y0(1)=BRAF_tot(k);
    y0(2)=ATP_tot(k);

    nexttile
    colororder(newcolors)

    [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), DBF_tot, TMT_tot, SCH_tot), tspan, y0, options);
    yo = y(:,38)/1.2;
    ypo = y(:,40)/1.2;
    yppo = y(:,44)/1.2;
    yo_pperk = y(:,39)/1.2;
    ypo_pperk = y(:,41)/1.2;
    yatp_o_pperk = y(:,42)/1.2;
    yatp_po_pperk = y(:,43)/1.2;
    ypo_phosph3 = y(:,46)/1.2;
    yppo_phosph3 = y(:,47)/1.2;
    yo_sch_pperk = y(:,50)/1.2;
    ypo_sch_pperk = y(:,51)/1.2;
    hold on
    plot(t,yo,t,ypo,t,yppo,t,yo_pperk,t,ypo_pperk, ...
        t,yatp_o_pperk,t,yatp_po_pperk,t,ypo_phosph3,t,yppo_phosph3, ...
        t,yo_sch_pperk,t,ypo_sch_pperk,'LineWidth',2.0)
    clear t y

    xticks(3600*[0 1 2 3])
    xlim([0 tspan(end)])
    xticklabels({'0','1','2','3'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)

    ylim([0 1])
    
    if k == 1
    xlabel('Time (hours)','FontSize',12)
    ylabel('Activity','FontSize',12)
    title('ATP = 1 mM', 'FontSize', 13)
    grid on
    grid minor

    elseif k == 2
    xlabel('Time (hours)','FontSize',12)
    ylabel('Activity','FontSize',12)
    annotation('textbox', [0.94, 0.65, 0.1, 0.1], 'String','BRAF = 3 nM', ...
    'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
    title('ATP = 5 mM', 'FontSize', 13)
    grid on
    grid minor

    elseif k == 3
    xlabel('Time (hours)','FontSize',12)
    ylabel('Activity','FontSize',12)
    grid on
    grid minor

    else
    xlabel('Time (hours)','FontSize',12)
    ylabel('Activity','FontSize',12)
    annotation('textbox', [0.945, 0.2, 0.1, 0.1], 'String','BRAF = 10 nM', ...
    'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
    grid on
    grid minor    
    end

end


lgd = legend('S', 'pS', 'ppS', 'S\cdot ppERK', 'pS\cdot ppERK', 'ATP\cdot S\cdot ppERK', 'ATP\cdot pS\cdot ppERK', ...
    'pS\cdot phosph3','ppS\cdot phosph3','S\cdot SCH\cdot ppERK','pS\cdot SCH\cdot ppERK', ...
    'Orientation', 'horizontal', 'Location','southoutside','Position',[0.5, 0.02, 1, 0.1]);
title(lgd, "Compound");
lgd.Layout.Tile = 'south';