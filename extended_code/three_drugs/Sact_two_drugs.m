close all
clear all


%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);
%Add a path to the directory that auxilary plotfiles
plot_location = strcat(pwd,'/auxiliary_files_plots');
addpath(plot_location);

DBF_min=0;
DBF_max=5;%10;
DBF_n=50;
DBF_h=(DBF_max-DBF_min)/DBF_n;
DBF_range=DBF_min:DBF_h:DBF_max;

TMT_min=0;
TMT_max=5;%10;
TMT_n=50;
TMT_h=(TMT_max-TMT_min)/TMT_n;
TMT_range=TMT_min:TMT_h:TMT_max;

SCH_min=0;
SCH_max=5;%2.5;
SCH_n=50;
SCH_h=(SCH_max-SCH_min)/SCH_n;
SCH_range=SCH_min:SCH_h:SCH_max;

timepoint_concentrations=zeros(length(DBF_range),length(DBF_range));
userinput_BRAF_tot = 0.003;
userinput_ATP_tot = 1000;
tend=[8 16 24];%in hours

set(gcf,'position',[800,1000,800,800])
tiles = tiledlayout(3,3,'TileSpacing','tight');
sp=1;
for row = 1:3

    if row == 1
        for col = 1:3
                for i_TMT = 1:length(TMT_range)
                    for i_DBF = 1:length(DBF_range)
                        timepoint_concentrations(i_TMT,i_DBF) = generate_concentration_data_sact(userinput_BRAF_tot, userinput_ATP_tot,...
                            DBF_range(i_DBF), TMT_range(i_TMT), SCH_min, tend(col));
                    end
                end
            nexttile, surf(DBF_range, TMT_range, timepoint_concentrations);       
            grid('on')
            xlabel('DBF (\muM)','FontSize',16)
            ylabel('TMT (\muM)' ,'FontSize',16)
            xlim([0 5])
            ylim([0 5])
            xticks([0 2 4])
            xticklabels({'0','2','4'})
            yticks([0 2 4])
            yticklabels({'0','2','4'})

            pbaspect([1 1 1]);
            view(2)
        
            if col==1
                title('8 h','FontSize',18)
            elseif col==2
                title('16 h','FontSize',18)
            elseif col==3
                title('24 h','FontSize',18)
            end
        end
    end

    if row == 2
        for col = 1:3
                for i_SCH = 1:length(SCH_range)
                    for i_DBF = 1:length(DBF_range)
                        timepoint_concentrations(i_SCH,i_DBF) = generate_concentration_data_sact(userinput_BRAF_tot, userinput_ATP_tot,...
                            DBF_range(i_DBF), TMT_min, SCH_range(i_SCH), tend(col));
                    end
                end
            nexttile, surf(DBF_range, SCH_range, timepoint_concentrations);       
            grid('on')
            xlabel('DBF (\muM)','FontSize',16)
            ylabel('SCH (\muM)' ,'FontSize',16)
            xlim([0 5])
            ylim([0 5])
            xticks([0 2 4])
            xticklabels({'0','2','4'})
            yticks([0 2 4])
            yticklabels({'0','2','4'})
            pbaspect([1 1 1]);
            view(2)
        
            if col==1
                title('8 h','FontSize',18)
            elseif col==2
                title('16 h','FontSize',18)
            elseif col==3
                title('24 h','FontSize',18)
            end
        end
    end

    if row == 3
        for col = 1:3
                for i_SCH = 1:length(SCH_range)
                    for i_TMT = 1:length(TMT_range)
                        timepoint_concentrations(i_SCH,i_TMT) = generate_concentration_data_sact(userinput_BRAF_tot, userinput_ATP_tot,...
                            DBF_min, TMT_range(i_TMT), SCH_range(i_SCH), tend(col));
                    end
                end
            nexttile, surf(TMT_range, SCH_range, timepoint_concentrations);       
            grid('on')
            xlabel('TMT (\muM)','FontSize',16)
            ylabel('SCH (\muM)' ,'FontSize',16)
            xlim([0 5])
            ylim([0 5])
            xticks([0 2 4])
            xticklabels({'0','2','4'})
            yticks([0 2 4])
            yticklabels({'0','2','4'})
            pbaspect([1 1 1]);
            view(2)
        
            if col==1
                title('8 h','FontSize',18)
            elseif col==2
                title('16 h','FontSize',18)
            elseif col==3
                title('24 h','FontSize',18)
            end
        end
    end

end
cb=colorbar('XTick',0:0.5:1, 'Location','eastoutside','Position',[0.93, 0.25, 0.02, 0.55]);
set(cb,'ylim',[0 1]);
cb.Title.String = {'Activated','substrate'};
cb.Title.FontSize = 14;