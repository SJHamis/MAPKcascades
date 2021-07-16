%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);
%Add a path to the directory that auxilary plotfiles
plot_location = strcat(pwd,'/auxiliary_files_plots');
addpath(plot_location);


DBF_min=0;
DBF_max=1;
DBF_n=10;
DBF_h=(DBF_max-DBF_min)/DBF_n;
DBF_range=DBF_min:DBF_h:DBF_max;

TMT_min=0;
TMT_max=1;
TMT_n=10;
TMT_h=(TMT_max-TMT_min)/TMT_n;
TMT_range=TMT_min:TMT_h:TMT_max;

%%%%%%%%%%%%%%%%%%%%%%%%% For investigating BRAF ranges
% figure
% threshold_braf=zeros(length(DBF_range),length(TMT_range));
% tend=24;
% 
% for i_TMT = 1:length(TMT_range)
%     for i_DBF = 1:length(DBF_range)
%         threshold_braf(i_TMT,i_DBF) = generate_concentration_data_brafv(...
%             DBF_range(i_DBF), TMT_range(i_TMT),tend);
%     end
% end
%         
% surf(DBF_range, TMT_range, threshold_braf);
% xlabel('DBF (\muM)','FontSize',14)
% ylabel('TMT (\muM)' ,'FontSize',14)
% zlabel('BRAF (nM)')
% 
% 
% set(gca,'ColorScale','log')
% cb=colorbar('eastoutside');
% cb.Ticks = [0.1 0.3 1 3 10 20]/1000;
% cb.TickLabels = {'0.1','0.3','1', '3', '10','20'};
% cb.Title.String ='BRAF (nM)';
% cb.Title.FontSize = 14;
% xticks([1 2 3 4 5])
% yticks([0.2 0.4 0.6 0.8 1])
% colormap(hot)
% view(0,90)
% cb.Title.FontSize = 14;
% 
% a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'fontsize',14)

% %%%%%%%%%%%%%%%%%%%%%%%%% For investigating ATP ranges

figure
threshold_atp=zeros(length(DBF_range),length(TMT_range));

tend=24;

for i_TMT = 1:length(TMT_range)
    for i_DBF = 1:length(DBF_range)
        threshold_atp(i_TMT,i_DBF) = generate_concentration_data_atpv(...
            DBF_range(i_DBF), TMT_range(i_TMT),tend);
    end
end
        
surf(DBF_range, TMT_range, threshold_atp);
xlabel('DBF (\muM)','FontSize',14)
ylabel('TMT (\muM)' ,'FontSize',14)
zlabel('ATP (mM)')

set(gca,'ColorScale','log')
cb=colorbar('eastoutside');
cb.Ticks = [0.1 0.3 1 3 10 20]*1000;
cb.TickLabels = {'0.1','0.3','1', '3', '10','20'};
cb.Title.String ='ATP (mM)';
cb.Title.FontSize = 14;
%xticks([1 2 3 4 5])
yticks([0 0.2 0.4 0.6 0.8 1])
xticks([0 0.2 0.4 0.6 0.8 1])
colormap(hot)
view(0,90)
cb.Title.FontSize = 14;

a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)