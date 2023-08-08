close all
clear all


%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);
%Add a path to the directory that auxilary plotfiles
plot_location = strcat(pwd,'/auxiliary_files_plots');
addpath(plot_location);

DBF_min=0;
DBF_max=5;
DBF_n=50;
DBF_h=(DBF_max-DBF_min)/DBF_n;
DBF_range=DBF_min:DBF_h:DBF_max;

TMT_min=0;
TMT_max=5;
TMT_n=50;
TMT_h=(TMT_max-TMT_min)/TMT_n;
TMT_range=TMT_min:TMT_h:TMT_max;

SCH_min=0;
SCH_max=5;
SCH_n=50;
SCH_h=(SCH_max-SCH_min)/SCH_n;
SCH_range=SCH_min:SCH_h:SCH_max;

SCH_concentration=zeros(length(TMT_range),length(DBF_range));

userinput_BRAF_tot = 0.003; % change to investigate other BRAF concentrations
userinput_ATP_tot = 1000; % change to investigate other ATP concentrations
tend=24; % in hours

set(gcf,'position',[100,100,1200,700])
    for i_SCH = 1:length(SCH_range)
            for i_DBF = 1:length(DBF_range)
                for i_TMT = 1:length(TMT_range)
                    timepoint_concentration = generate_concentration_data_sact(userinput_BRAF_tot, userinput_ATP_tot,...
                        DBF_range(i_DBF), TMT_range(i_TMT), SCH_range(i_SCH), tend);
                    if i_SCH == 1 && timepoint_concentration<=0.50  %change for IC25, IC50 or IC75
                            SCH_concentration(i_DBF, i_TMT) = 100; %dummy value to preserve data for SCH = 0
                    elseif timepoint_concentration<=0.50 && SCH_concentration(i_DBF, i_TMT) == 0
                        SCH_concentration(i_DBF, i_TMT) = SCH_range(i_SCH);
                    elseif i_SCH == 51 && SCH_concentration(i_DBF,i_TMT) == 100
                        SCH_concentration(i_DBF,i_TMT) = 0; %get SCH=0 back at the end of the loop
                    end
                end
            end
    end

surf(TMT_range, DBF_range, SCH_concentration);
a = get(gca,'Clim');

hold on

% plane eqn through 3 monotherapy points calculated separately

% A = 2.89 for BRAF = 3 nM, ATP = 1 mM;
%   = 10.2 for BRAF = 10 nM, ATP = 1 mM;
%   = 37.45 for BRAF = 3 nM, ATP = 5 mM.
A = 2.89;

% B = 4.93 for BRAF = 3 nM, ATP = 1 mM;
%   = 4.93 for BRAF = 10 nM, ATP = 5 mM;
%   = 18.2 for BRAF = 3 nM, ATP = 1 mM.
B = 4.93;

% C = 4.93 for BRAF = 3 nM, ATP = 1 mM;
%   = 17.4 for BRAF = 10 nM, ATP = 1 mM;
%   = 55.64 for BRAF = 3 nM, ATP = 5 mM.
C = 4.93;

% D = -8.381 for BRAF = 3 nM, ATP = 1 mM;
%   = -29.58 for BRAF = 10 nM, ATP = 1 mM;
%   = -194.74 for BRAF = 3 nM, ATP = 5 mM.
D = -8.381;

[X,Y] = meshgrid(TMT_range,DBF_range);
synergy_plane = (-1/C)*(A*X + B*Y + D); % synergy plane
surf(TMT_range, DBF_range, synergy_plane, ...
    'FaceColor','red', 'FaceAlpha',0.5, 'EdgeColor','none');

axis([0 5 0 5 0 5])
set(gca,'Clim',a);

grid on
xlabel('TMT (\muM)','FontSize',14)
ylabel('DBF (\muM)' ,'FontSize',14)
set(gca, 'YDir','reverse')
zlabel('SCH (\muM)', 'FontSize',14)

title({'IC_{50} : BRAF = 3 nM, ATP = 1 mM'}, 'FontSize', 14)