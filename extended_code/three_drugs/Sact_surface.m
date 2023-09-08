close all
clear all


%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);
%Add a path to the directory that auxilary plotfiles
plot_location = strcat(pwd,'/auxiliary_files_plots');
addpath(plot_location);

DBF_min=0;
DBF_max=2;
DBF_n=40;
DBF_h=(DBF_max-DBF_min)/DBF_n;
DBF_range=DBF_min:DBF_h:DBF_max;

TMT_min=0;
TMT_max=2;
TMT_n=40;
TMT_h=(TMT_max-TMT_min)/TMT_n;
TMT_range=TMT_min:TMT_h:TMT_max;

SCH_min=0;
SCH_max=2;
SCH_n=40;
SCH_h=(SCH_max-SCH_min)/SCH_n;
SCH_range=SCH_min:SCH_h:SCH_max;

SCH_concentration=zeros(length(TMT_range),length(DBF_range));

userinput_BRAF_tot = 0.003; % change to investigate other BRAF concentrations
userinput_ATP_tot = 1000; % change to investigate other ATP concentrations
tend = 8; % in hours
IC_in = 25; % change for IC25, IC50, and IC75


% draw three-drug combination surface

set(gcf,'position',[100,100,1200,700])
    for i_SCH = 1:length(SCH_range)
            for i_DBF = 1:length(DBF_range)
                for i_TMT = 1:length(TMT_range)
                    timepoint_concentration = generate_concentration_data_sact(userinput_BRAF_tot, userinput_ATP_tot,...
                        DBF_range(i_DBF), TMT_range(i_TMT), SCH_range(i_SCH), tend);
                    if i_SCH == 1 && timepoint_concentration <= IC_in/100
                            SCH_concentration(i_DBF, i_TMT) = 100; %dummy value to preserve data for SCH = 0
                    elseif timepoint_concentration <= IC_in/100 && SCH_concentration(i_DBF, i_TMT) == 0
                        SCH_concentration(i_DBF, i_TMT) = SCH_range(i_SCH);
                    elseif i_SCH == (SCH_n+1) && SCH_concentration(i_DBF,i_TMT) == 100
                        SCH_concentration(i_DBF,i_TMT) = 0; %get SCH=0 back at the end of the loop
                    end
                end
            end
    end

surf(TMT_range, DBF_range, SCH_concentration);
ax = get(gca,'Clim');

hold on


% add synergy plane

% find monotherapy values for the defined IC and for each drug
monotherapy_points = generate_monotherapy_dose_data(userinput_BRAF_tot, userinput_ATP_tot, IC_in, tend);

% choose one point lying on the plane P = (x0, y0, z0)
P = monotherapy_points(1,:);
x0 = P(1);
y0 = P(2);
z0 = P(3);

% compute the normal vector n = (a, b, c)
i = monotherapy_points(2,:)-P;
j = monotherapy_points(3,:)-P;
n = cross(i,j);
a = n(1);
b = n(2);
c = n(3);

% use the formula to find the plane passing through the three monotherapy points
% a(x-x0)+b(y-y0)+c(z-z0)=0 -> A = a; B = b; C = c; D = -a*x0
A = a;
B = b;
C = c;
D = -a*x0;

[X,Y] = meshgrid(TMT_range,DBF_range);
synergy_plane = (-1/C)*(A*X + B*Y + D); % synergy plane
surf(TMT_range, DBF_range, synergy_plane, ...
    'FaceColor','red', 'FaceAlpha',0.5, 'EdgeColor','none');

axis([0 2 0 2 0 2])
set(gca,'Clim',ax);

grid on
xlabel('TMT (\muM)','FontSize',14)
ylabel('DBF (\muM)' ,'FontSize',14)
set(gca, 'YDir','reverse')
zlabel('SCH (\muM)', 'FontSize',14)

title({'IC_{25} : BRAF = 3 nM, ATP = 5 mM; 8 h'}, 'FontSize', 14)