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
MEK_tot=1.2;
ERK_tot=1.2;
SUB_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
phosph3_tot=0.12;

y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot;
y0(38)=SUB_tot;
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;
y0(42)=phosph3_tot;

DBF_min = 0;
DBF_max = 20;
n_DBF = 1000;
DBF_vector = linspace(DBF_min, DBF_max, n_DBF);

TMT_min = 0;
TMT_max = 20;
n_TMT = 1000;
TMT_vector = linspace(TMT_min, TMT_max, n_TMT);

SCH_min = 0;
SCH_max = 20;
n_SCH = 1000;
SCH_vector = linspace(SCH_min, SCH_max, n_SCH);

options1 = odeset('Mass',M, 'MassSingular','yes', 'RelTol',1e-4, 'AbsTol',1e-4);
options2 = odeset('Mass',M, 'MassSingular','yes', 'RelTol',3e-4, 'AbsTol',3e-4);
set(gcf,'Position',[100,100,1000,1000])
newcolors = {'#EDB120','#7E2F8E','#77AC30','#D95319','#0072BD','#FF00FF'};

%%%%%%%%%%%%%%%%%%%%%%% Generate activated substrate SUB (unbound pSUB) vs monotherapy (DBF, TMT, SCH) FIGURE for different t values
tiles = tiledlayout(3,4, TileSpacing="compact", Padding="loose"); 

BRAF_tot=[0.003,0.01,0.003,0.01, 0.003,0.01,0.003,0.01, 0.003,0.01,0.003,0.01];
ATP_tot=[1000,1000,5000,5000, 1000,1000,5000,5000, 1000,1000,5000,5000];

for k = 1:12 
    y0(1)=BRAF_tot(k);
    y0(2)=ATP_tot(k);

    nexttile
    colororder(newcolors)


    if k <= 4  % plotting for DBF monotherapy

        output_data = zeros(6,length(DBF_vector));   % storage matrix to plot results
        y0(30)=0;
        y0(44)=0;
        index = 1;  % to update the output_data matrix

        for i = [3 8 16 24 60 72]
            tspan = [0, i*3600];

            for j = 1:length(DBF_vector)

                y0(15)=DBF_vector(j);
                DBF_in = y0(15);
            
                [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), DBF_in, 0, 0), tspan, y0, options1);
                ysubact = y(:,41)/1.2;

                output_data(index,j) = ysubact(end);
            end

            index = index+1;
        end

        plot(DBF_vector,output_data(1,:),'LineWidth',6);  % plot after 3 hours
        hold on
        plot(DBF_vector,output_data(2,:),'LineWidth',5);  % plot after 8 hours
        hold on
        plot(DBF_vector,output_data(3,:),'LineWidth',4);  % plot after 16 hours
        hold on
        plot(DBF_vector,output_data(4,:),'LineWidth',3);  % plot after 24 hours
        hold on
        plot(DBF_vector,output_data(5,:),'LineWidth',2);  % plot after 60 hours
        hold on
        plot(DBF_vector,output_data(6,:),'LineWidth',1);  % plot after 72 hours
        clear t y
    
        xticks([0 5 10 15 20])
        xlim([0 DBF_vector(end)])
        xticklabels({'0','5','10','15','20'})
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12)
        pbaspect([1 1 1])
    
        ylim([0 0.3])
    
        xlabel('DBF dose (\muM)','FontSize',12)
        ylabel('Activated substrate','FontSize',12)
        grid on
        grid minor
        
        if k == 1
            title({'BRAFV600E_{tot} = 3 nM,', 'ATP_{tot} = 1 mM'}, 'FontSize', 12)
            annotation('textbox', [0.05, 0.76, 0.1, 0.1], 'String','DBF', ...
            'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);    
    
        elseif k == 2
            title({'BRAFV600E_{tot} = 10 nM,', 'ATP_{tot} = 1 mM'}, 'FontSize', 12)
    
        elseif k == 3
            title({'BRAFV600E_{tot} = 3 nM,', 'ATP_{tot} = 5 mM'}, 'FontSize', 12)
    
        elseif k == 4
            title({'BRAFV600E_{tot} = 10 nM,', 'ATP_{tot} = 5 mM'}, 'FontSize', 12)

        end

    end


    if k > 4 && k <=8 % plotting for TMT monotherapy

        output_data = zeros(6,length(TMT_vector));   % storage matrix for plot results
        y0(15)=0;
        y0(44)=0;
        index = 1;  % to update the output_data matrix

        for i = [3 8 16 24 60 72]
            tspan = [0, i*3600];

            for j = 1:length(TMT_vector)

                y0(30)=TMT_vector(j);
                TMT_in = y0(30);
            
                [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), 0, TMT_in, 0), tspan, y0, options2);
                ysubact = y(:,41)/1.2;

                output_data(index,j) = ysubact(end);
            end

            index = index+1;
        end

        plot(TMT_vector,output_data(1,:),'LineWidth',6);
        hold on
        plot(TMT_vector,output_data(2,:),'LineWidth',5);
        hold on
        plot(TMT_vector,output_data(3,:),'LineWidth',4);
        hold on
        plot(TMT_vector,output_data(4,:),'LineWidth',3);
        hold on
        plot(TMT_vector,output_data(5,:),'LineWidth',2);
        hold on
        plot(TMT_vector,output_data(6,:),'LineWidth',1);
        clear t y
    
        xticks([0 5 10 15 20])
        xlim([0 TMT_vector(end)])
        xticklabels({'0','5','10','15','20'})
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12)
        pbaspect([1 1 1])
    
        ylim([0 0.3])
    
        xlabel('TMT dose (\muM)','FontSize',12)
        ylabel('Activated substrate','FontSize',12)
        grid on
        grid minor
        
        if k == 5
            annotation('textbox', [0.05, 0.47, 0.1, 0.1], 'String','TMT', ...
            'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);    

        end

    end


    if k > 8 && k <=12 % plotting for SCH772984 monotherapy

        output_data = zeros(6,length(SCH_vector));   % storage matrix for plot results
        y0(15)=0;
        y0(30)=0;
        index = 1;  % to update the output_data matrix

        for i = [3 8 16 24 60 72]
            tspan = [0, i*3600];

            for j = 1:length(SCH_vector)

                y0(44)=SCH_vector(j);
                SCH_in = y0(44);
            
                [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), 0, 0, SCH_in), tspan, y0, options1);
                ysubact = y(:,41)/1.2;
       
                output_data(index,j) = ysubact(end);
            end

            index = index+1;
        end
        
        hold on
        plot(SCH_vector,output_data(1,:),'LineWidth',6);
        hold on
        plot(SCH_vector,output_data(2,:),'LineWidth',5);
        hold on
        plot(SCH_vector,output_data(3,:),'LineWidth',4);
        hold on
        plot(SCH_vector,output_data(4,:),'LineWidth',3);
        hold on
        plot(SCH_vector,output_data(5,:),'LineWidth',2);
        hold on
        plot(SCH_vector,output_data(6,:),'LineWidth',1);
        clear t y
    
        xticks([0 5 10 15 20])
        xlim([0 SCH_vector(end)])
        xticklabels({'0','5','10','15','20'})
        a = get(gca,'XTickLabel');
        set(gca,'XTickLabel',a,'fontsize',12)
        pbaspect([1 1 1])
    
        ylim([0 0.3])
    
        xlabel('SCH dose (\muM)','FontSize',12)
        ylabel('Activated substrate','FontSize',12)
        grid on
        grid minor
        
        if k == 9
            annotation('textbox', [0.05, 0.18, 0.1, 0.1], 'String','SCH', ...
            'EdgeColor','none', 'FontWeight', 'bold', 'FontSize', 14, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);    

        end

    end

end

lgd = legend('3 h','8 h','16 h','24 h','60 h','72 h', ...
    'Orientation', 'horizontal', 'Location','south');
title(lgd, "Time (hours)");
lgd.Layout.Tile = 'south';