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

y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot;
y0(38)=O_tot;
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;
y0(45)=phosph3_tot;
y0(15)=DBF_tot;
y0(30)=TMT_tot;

SCH_min = 0;
SCH_max = 10;
n_SCH = 1000;
SCH_vector = linspace(SCH_min, SCH_max, n_SCH);

options = odeset('Mass',M, 'MassSingular','yes', 'RelTol',1e-3, 'AbsTol',1e-3);
set(gcf,'position',[100,100,1200,700])
newcolors = {'#FF5733 ','#0072BD','#00FF00','#008000',' #FFC300 ','#900C3F','#581845','#000000','#DAF7A6'};

%%%%%%%%%%%%%%%%%%%%%%% Generate activated substrate O (pO + ppO) vs SCH772984 FIGURE for different t values
tiles = tiledlayout(2,2,'TileSpacing','compact'); 

BRAF_tot=[0.003,0.003,0.01,0.01];
ATP_tot=[1000, 5000, 1000, 5000];

for k = 1:4 
    y0(1)=BRAF_tot(k);
    y0(2)=ATP_tot(k);

    nexttile
    colororder(newcolors)

    output_data = zeros(3,length(SCH_vector));   % storage matrix for plot results
    
    for i = [1 2 3]
        tspan = [0, i*3600];

        for j = 1:length(SCH_vector)
            y0(48)=SCH_vector(j);
            SCH_in = y0(48);
            
            [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), DBF_tot, TMT_tot, SCH_in), tspan, y0, options);
            ypo = y(:,40)/1.2;
            yppo = y(:,44)/1.2;
            yoact = (ypo+yppo);

            output_data(i,j) = yoact(end);
        end
    end
    
    hold on
    plot(SCH_vector,output_data(1,:),SCH_vector,output_data(2,:),SCH_vector,output_data(3,:), ...
        'LineWidth',2.0)
    clear t y

    xticks([0 1 2 3 4 5 6 7 8 9 10])
    xlim([0 SCH_vector(end)])
    xticklabels({'0','1','2','3','4','5','6','7','8','9','10'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',12)

    ylim([0 1])
    
    if k == 1
    xlabel('SCH dose (\muM)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    title('ATP = 1 mM', 'FontSize', 13)
    grid on
    grid minor

    elseif k == 2
    xlabel('SCH dose (\muM)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    annotation('textbox', [0.94, 0.65, 0.1, 0.1], 'String','BRAF = 3 nM', ...
    'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
    title('ATP = 5 mM', 'FontSize', 13)
    grid on
    grid minor

    elseif k == 3
    xlabel('SCH dose (\muM)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    grid on
    grid minor

    else
    xlabel('SCH dose (\muM)','FontSize',12)
    ylabel('Predicted activated substrate activity','FontSize',12)
    annotation('textbox', [0.945, 0.2, 0.1, 0.1], 'String','BRAF = 10 nM', ...
    'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
    grid on
    grid minor    
    end

end


lgd = legend('1 h','2 h','3 h', ...
    'Orientation', 'horizontal', 'Location','southoutside','Position',[0.5, 0.02, 1, 0.1]);
title(lgd, "Time (hours)");
lgd.Layout.Tile = 'south';