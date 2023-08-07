close all
clear all


%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);
%Add a path to the directory that auxilary plotfiles
plot_location = strcat(pwd,'/auxiliary_files_plots');
addpath(plot_location);

DBF_in=0;

TMT_min=0;
TMT_max=5;%10;
TMT_n=25;
TMT_h=(TMT_max-TMT_min)/TMT_n;
TMT_range=TMT_min:TMT_h:TMT_max;


SCH_min=0;
SCH_max=5;%2.5;
SCH_n=25;
SCH_h=(SCH_max-SCH_min)/SCH_n;
SCH_range=SCH_min:SCH_h:SCH_max;

timepoint_concentrations=zeros(length(TMT_range),length(SCH_range));
userinput_BRAF_tot =[0.003, 0.01, 0.003];
userinput_ATP_tot = [1000, 1000, 5000];
tend=[8 16 24];%in hours

set(gcf,'position',[100,100,1200,700])
tiles = tiledlayout(3,3,'TileSpacing','compact');
sp=1;
for row = 1:3
    for col = 1:3
            for i_SCH = 1:length(SCH_range)
                for i_TMT = 1:length(TMT_range)
                    timepoint_concentrations(i_SCH,i_TMT) = generate_concentration_data_sact(userinput_BRAF_tot(row), userinput_ATP_tot(row),...
                        DBF_in, TMT_range(i_TMT), SCH_range(i_SCH), tend(col));
                end
            end
        nexttile, surf(TMT_range, SCH_range, timepoint_concentrations);       
        grid('on')
        xlabel('TMT (\muM)','FontSize',14)
        ylabel('SCH (\muM)' ,'FontSize',14)      
        view(2)

        if row==1 && col==1
            title('8 h')
            annotation('textbox', [0.089, 0.65, 0.1, 0.1], 'String','BRAF = 3 nM, ATP = 1 mM', ...
            'EdgeColor','none', 'FontSize', 10, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
        elseif row==1 && col==2
            title('16 h')
        elseif row==1 && col==3
            title('24 h')
        elseif row==2 && col==1 
            annotation('textbox', [0.09, 0.35, 0.1, 0.1], 'String','BRAF = 10 nM, ATP = 1 mM', ...
            'EdgeColor','none', 'FontSize', 10, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
        elseif row==3 && col==1
            annotation('textbox', [0.088, 0.05, 0.1, 0.1], 'String','BRAF = 3 nM, ATP = 5 mM', ...
            'EdgeColor','none', 'FontSize', 10, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment','middle', 'Rotation',90);
        end

    end
end
cb=colorbar('XTick',0:0.5:1, 'Location','eastoutside','Position',[0.93, 0.25, 0.02, 0.5]);
set(cb,'ylim',[0 1]);
cb.Title.String = {'Activated','substrate'}
cb.Title.FontSize = 12;