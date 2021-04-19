close all
clear all

%Add a path to the directory that contains the model details
dae_location = strcat(pwd,'/auxiliary_files_model_setup');
addpath(dae_location);

M = eye(36);
super_compound_list_index=get_conslaw_position();
%Substitute in cons. laws 
for i=1:size(M,1)
    if(ismember(i,super_compound_list_index))
        M(i,i)=0; 
    end
end

% Set initial conditions
MEK_tot=1.2;
ERK_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot; 
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;

DBFlist=[0 1 5 10 50 100]/10;
TMTlist=DBFlist;%/10;
tend = [8 16 24]*3600;

options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-4,'AbsTol',1e-4);
set(gcf,'position',[100,100,1000,600])
newcolors = {'#008000','#008000','#008000',...
    '#FFC300 ','#FFC300 ','#FFC300 ',...
    '#FF5733 ','#FF5733 ','#FF5733 ',...
    '#C70039','#C70039','#C70039',...
    '#581845','#581845','#581845',...
    '#000000','#000000','#000000'};
    
%%%%%%%%%%%%%%%%%%%%%%% FIG 1
t = tiledlayout(2,2,'TileSpacing','none','Padding','none');
TMT_tot=0;
y0(30)=TMT_tot;
TMT_in=TMT_tot;

%%%this plot
BRAF_range=linspace(0.003/10,0.01,100);
ATP_tot=1000;
y0(2)=ATP_tot;
output = zeros(length(BRAF_range), length(DBFlist), length(tend));
%%%%%%


%loop through all the braf values
for k_braf = 1:length(BRAF_range)
    y0(1)=BRAF_range(k_braf);   
    BRAF_tot = BRAF_range(k_braf);   
        
    %loop through all dbf concentrations
    for j_dbf = 1:length(DBFlist)
        DBF_in=DBFlist(j_dbf);
        y0(15)=DBF_in;
        %output(k_braf, 2) = DBF_in;
        
        %loop through all 3 end timepoints
        for l_endt = 1:length(tend)
            tspan = [0 tend(l_endt)];
            
            [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot, ATP_tot, DBF_in, TMT_in), tspan, y0, options);
            perkpperk = y(end,26)/1.2;                
            output(k_braf,j_dbf,l_endt) = perkpperk; 
        end
    end
end
 
nexttile
colororder(newcolors)
for j_dbf=1:length(DBFlist)
    hold on
    plot(BRAF_range, output(:,j_dbf,1),'--', BRAF_range, output(:,j_dbf,2),':', BRAF_range, output(:,j_dbf,3),'LineWidth',2)
end
   
ylabel('Activated ERK','FontSize',18)
yticks([0 0.5 1])  

xlabel('BRAF (nM)','FontSize',18)
xticks([BRAF_range(1) 0.003 0.010])    
xticklabels({'0.3','3','10'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)


xlim([BRAF_range(1), BRAF_range(end)])
grid on
    
%%%%% ATP range 
ATP_range=linspace(1000,5000,100);
BRAF_tot=0.003;
y0(1)=BRAF_tot;
output = zeros(length(ATP_range), length(DBFlist), length(tend));

%loop through all the atp values
for k_atp = 1:length(ATP_range)
    y0(2)=ATP_range(k_atp);   
    ATP_tot = ATP_range(k_atp);   
        
    %loop through all dbf concentrations
    for j_dbf = 1:length(DBFlist)
        DBF_in=DBFlist(j_dbf);
        y0(15)=DBF_in;      
        
        %loop through all 3 end timepoints
        for l_endt = 1:length(tend)
            tspan = [0 tend(l_endt)];
            
            [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot, ATP_tot, DBF_in, TMT_in), tspan, y0, options);
            perkpperk = y(end,26)/1.2;                
            output(k_atp,j_dbf,l_endt) = perkpperk; 
        end
    end
end

nexttile
colororder(newcolors)
for j_dbf=1:length(DBFlist)
    hold on
    plot(ATP_range, output(:,j_dbf,1),'--', ATP_range, output(:,j_dbf,2),':', ATP_range, output(:,j_dbf,3),'LineWidth',2)
end

ylabel('Activated ERK','FontSize',18)
yticks([0 0.5 1])
xlabel('ATP (mM)','FontSize',18)
xlim([ATP_range(1), ATP_range(end)])
xticks([1000 3000 5000])    
xticklabels({'1','3','5'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
grid on   

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FIG 2 TMT
options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-3,'AbsTol',1e-3);
DBF_tot=0;
y0(15)=DBF_tot;
DBF_in=DBF_tot;

%%%this plot
BRAF_range=linspace(0.003/10,0.010,100);
ATP_tot=1000;
y0(2)=ATP_tot;
output = zeros(length(BRAF_range), length(TMTlist), length(tend));
%%%%%%


%loop through all the braf values
for k_braf = 1:length(BRAF_range)
    y0(1)=BRAF_range(k_braf);   
    BRAF_tot = BRAF_range(k_braf);   
        
    %loop through all tmt concentrations
    for j_tmt = 1:length(TMTlist)
        TMT_in=TMTlist(j_tmt);
        y0(30)=TMT_in;
        %output(k_braf, 2) = DBF_in;
        
        %loop through all 3 end timepoints
        for l_endt = 1:length(tend)
            tspan = [0 tend(l_endt)];
            
            [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot, ATP_tot, DBF_in, TMT_in), tspan, y0, options);
            perkpperk = (y(end,26))/1.2;                
            output(k_braf,j_tmt,l_endt) = perkpperk; 
        end
    end
end
 
nexttile
colororder(newcolors)
for j_tmt=1:length(TMTlist)
    hold on
    plot(BRAF_range, output(:,j_tmt,1),'--', BRAF_range, output(:,j_tmt,2),':', BRAF_range, output(:,j_tmt,3),'LineWidth',2)
end
ylabel('Activated ERK','FontSize',18)
yticks([0 0.5 1])  

xlabel('BRAF (nM)','FontSize',18)
xlim([BRAF_range(1), BRAF_range(end)])

xticks([BRAF_range(1) 0.003 0.01])
xticklabels({'0.3','3','10'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
grid on
    

ATP_range=linspace(1000,5000,100);
BRAF_tot=0.003;
y0(1)=BRAF_tot;
output = zeros(length(ATP_range), length(TMTlist), length(tend));

%loop through all the braf values
for k_atp = 1:length(ATP_range)
    y0(2)=ATP_range(k_atp);   
    ATP_tot = ATP_range(k_atp);   
        
    %loop through all tmt concentrations
    for j_tmt = 1:length(TMTlist)
        TMT_in=TMTlist(j_tmt);
        y0(30)=TMT_in;      
        
        %loop through all 3 end timepoints
        for l_endt = 1:length(tend)
            tspan = [0 tend(l_endt)];
            
            [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot, ATP_tot, DBF_in, TMT_in), tspan, y0, options);
            perkpperk = (y(end,26))/1.2;                
            output(k_atp,j_tmt,l_endt) = perkpperk; 
        end
    end
end
 
nexttile
colororder(newcolors)
for j_tmt=1:length(TMTlist)
    hold on
    plot(ATP_range, output(:,j_tmt,1),'--', ATP_range, output(:,j_tmt,2),':', ATP_range, output(:,j_tmt,3),'LineWidth',2)
end

ylabel('Activated ERK','FontSize',18)
yticks([0 0.5 1])
xlabel('ATP (mM)','FontSize',18)
xlim([ATP_range(1), ATP_range(end)])
xticks([1000 3000 5000])    
xticklabels({'1','3','5'})
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)
   
grid on
  
    
