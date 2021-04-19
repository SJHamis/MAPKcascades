%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);
%Add a path to the directory that auxilary plotfiles
plot_location = strcat(pwd,'/auxilary_files_plots');
addpath(plot_location);

DBF_min=0;
DBF_max=2;%10;
DBF_n=20;
DBF_h=(DBF_max-DBF_min)/DBF_n;
DBF_range=DBF_min:DBF_h:DBF_max;

TMT_min=0;
TMT_max=2;%2.5;
TMT_n=20;
TMT_h=(TMT_max-TMT_min)/TMT_n;
TMT_range=TMT_min:TMT_h:TMT_max;

timepoint_concentrations=zeros(length(DBF_range),length(TMT_range));
userinput_BRAF_tot =[0.003, 0.01, 0.003];
userinput_ATP_tot = [1000, 1000, 5000];
tend=[8 16 24];%in hours

t = tiledlayout(3,3,'TileSpacing','none','Padding','none');
sp=1;
for row = 1:3
    for col = 1:3       
        for i_TMT = 1:length(TMT_range)
            for i_DBF = 1:length(DBF_range)
                timepoint_concentrations(i_TMT,i_DBF) = generate_concentration_data(userinput_BRAF_tot(row), userinput_ATP_tot(row),...
                    DBF_range(i_DBF), TMT_range(i_TMT),tend(col));
            end
        end
        nexttile, surf(DBF_range, TMT_range, timepoint_concentrations);       
        grid('on')
        xlabel('DBF (\muM)','FontSize',14)
        ylabel('TMT (\muM)' ,'FontSize',14)      
        view(2)
    end
end

cb=colorbar('southoutside', 'XTick',0:0.5:1);
set(cb,'ylim',[0 1]);
cb.Title.String = 'Activated ERK';
cb.Title.FontSize = 14;
     
