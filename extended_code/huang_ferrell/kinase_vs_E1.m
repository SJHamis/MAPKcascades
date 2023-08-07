close all
clear all

%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);

% Set inital conditions
KKK_tot=0.003;
KK_tot=1.2;
K_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
E2_tot=0.0003;
y0=zeros(22,1);
y0(1)=KKK_tot;
y0(5)=E2_tot;
y0(7)=KK_tot;
y0(15)=K_tot; 
y0(12)=phosph1_tot;
y0(20)=phosph2_tot;

KKK_in=KKK_tot
E2_in=E2_tot

E1_min = 0.000001;
E1_max = 0.1;
n_E1 = 200000;
E1_vector = linspace(E1_min, E1_max, n_E1);
t_span = [456*3600, 480*3600]

set(gcf,'position',[100,100,1000,600])
newcolors = {'#008000',...
    '#FFC300',...
    '#FF5733 ',...
    '#C70039','#C70039','#C70039',...
    '#581845','#581845','#581845',...
    '#000000','#000000','#000000'};
    
%%%%%%%%%%%%%%%%%%%%%%% Generate figure activated kinases vs E1

output_data = zeros(3,length(E1_vector));   %storage matrix for plot results

%loop through all the E1 values
for j = 1:length(E1_vector)
    y0(2)=E1_vector(j);   
    E1_in = y0(2);

    [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, KKK_in, E1_in, E2_in), t_span, y0);
    y_plot_KKK = y(:,4) ./ (y(:,1)+y(:,4));
    y_plot_KK = y(:,11) ./ (y(:,7)+y(:,9)+y(:,11));
    y_plot_K = y(:,19) ./ (y(:,15)+y(:,17)+y(:,19));

    output_data(1,j) = y_plot_KKK(end);
    output_data(2,j) = y_plot_KK(end);
    output_data(3,j) = y_plot_K(end);
end
 
colororder(newcolors)
semilogx(E1_vector, output_data(1,:), E1_vector, output_data(2,:), E1_vector, output_data(3,:),'LineWidth',1.0)
yline(0.1, '--', "10%" )
yline(0.9, "--", "90%")

ylabel('Predicted kinase activity','FontSize',12) 
xlabel('E1 (nM)','FontSize',12)
title('Activated Kinase vs. E1', 'FontSize', 13)
legend("KKK*", "KK-PP", "K-PP", ...
    "Location", "east")

grid on
grid minor