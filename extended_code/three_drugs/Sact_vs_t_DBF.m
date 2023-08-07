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
TMT_tot=0;
SCH_tot=0;

y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot;
y0(38)=O_tot;
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;
y0(45)=phosph3_tot;
y0(30)=TMT_tot;
y0(48)=SCH_tot;


tspan = [0 48*3600];   % 48h (can be changed later)
DBF_tot=[0 1 5 10 25 50 100]/10;

options = odeset('Mass',M, 'MassSingular','yes', 'RelTol',5e-4, 'AbsTol',5e-4);
set(gcf,'position',[100,100,1200,700])
newcolors = {'#008000',' #FFC300 ', '#FF5733 ','#C70039','#900C3F','#581845','#000000','#DAF7A6'};

%%%%%%%%%%%%%%%%%%%%%%% Generate activated O (pO + ppO) vs t FIGURE for different DBF values
tiles = tiledlayout(2,2,'TileSpacing','compact'); 

BRAF_tot=[0.003,0.003,0.01,0.01];
ATP_tot=[1000, 5000, 1000, 5000];

for k = 1:4
    y0(1)=BRAF_tot(k);
    y0(2)=ATP_tot(k);

    nexttile
    colororder(newcolors)

    for j = 1:length(DBF_tot)
        DBF_in=DBF_tot(j);
        y0(15)=DBF_in;
        [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), DBF_in, TMT_tot, SCH_tot), tspan, y0, options);
        ypo = y(:,40)/1.2;
        yppo = y(:,44)/1.2;
        yoact = ypo+yppo;
        hold on
        plot(t,yoact,'LineWidth',2.0)
        clear t y
    end

    xticks(3600*[0 8 16 24 48])
    xlim([0 tspan(end)])
    xticklabels({'0','8','16','24','48'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)

    ylim([0 1])
    
    if k == 1
    xlabel('Time (hours)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    title('ATP = 1 mM', 'FontSize', 13)
    grid on
    grid minor

    elseif k == 2
    xlabel('Time (hours)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    annotation('textbox', [0.94, 0.65, 0.1, 0.1], 'String','BRAF = 3 nM', ...
    'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
    title('ATP = 5 mM', 'FontSize', 13)
    grid on
    grid minor

    elseif k == 3
    xlabel('Time (hours)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    grid on
    grid minor

    else
    xlabel('Time (hours)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    annotation('textbox', [0.945, 0.2, 0.1, 0.1], 'String','BRAF = 10 nM', ...
    'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
    grid on
    grid minor    
    end


end

lgd = legend('0\muM', '0.1\muM', '0.5\muM', '1\muM', '2.5\muM', '5\muM', '10\muM', ...
'Orientation', 'horizontal', 'Location','southoutside','Position',[0.5, 0.02, 1, 0.1]);
title(lgd, "DBF doses");
lgd.Layout.Tile = 'south';