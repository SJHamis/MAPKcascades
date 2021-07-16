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

% Set non-changing inital conditions
MEK_tot=1.2;
ERK_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot; 
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;
TMT_tot=0;
y0(30)=TMT_tot;
TMT_in=TMT_tot;


dbfdose = linspace(0,1.9246,100);
dbfresponse=zeros(1,length(dbfdose));

options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-4,'AbsTol',1e-4);
set(gcf,'position',[100,100,1000,600])
newcolors = {'#008000',' #FFC300 ', '#FF5733 ','#C70039','#900C3F','#581845','#000000','#DAF7A6'};


colororder(newcolors)

%Input: 1:Hour, 2:BRAF factor aplifiction, 3:ATP factor amplification  
inputM=[8, 1, 1, 0;
        8, 3.3333, 1, 0];    
   
setEndTime=3600*inputM(:,1);
setBRAF=0.003*inputM(:,2);
setATP=1000*inputM(:,3);

for v = 1:2%length(inputM)
    tspan = [0 setEndTime(v)];    
    BRAF_tot=setBRAF(v);
    y0(1)=BRAF_tot;
    ATP_tot=setATP(v);        
    y0(2)=ATP_tot;

    for d = 1:length(dbfdose)
        DBF_in=dbfdose(d);
        y0(15)=DBF_in; 
        [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot, ATP_tot, DBF_in, TMT_in), tspan, y0, options);
        yperk = y(end,22);
        ypperk = y(end,26);
        dbfresponse(d)=(yperk+ypperk)/1.2;
    end

    dbfresponse=dbfresponse/dbfresponse(1);
    hold on
    if(v==1)
        plot(dbfdose,dbfresponse,'k','LineWidth',2)
    else
        plot(dbfdose,dbfresponse,'k--','LineWidth',2)
    end
end

% uMol -> ng conversion:
%[ng]=[mol/liter]*[mL]*[g/mol]*1000(to make it nano)
%effconc = dbfdose*0.001*519.6*1000;
%[0:200:1000]/(0.001*519.6*1000);


yticks([0:0.2:1])
yticklabels({'-100','-80','-60','-40','-20','0'})
ylim([0 1])
ylabel('Change in phosphorylated ERK (%)','FontSize',18)

% XAxisLocation = 'bottom';
xticks([0:200:1000]/(0.001*519.6*1000))
xticklabels({'0','200','400','600','800','1000'})
xlim([0 dbfdose(end)])
xlabel('Effective dabrafenib concentration (ng/mL)','FontSize',18)


%Fetching data from Falchook et al.'s plot via:
%https://apps.automeris.io/wpd/
data=[137.3913043478261, 62.32365145228215
177.97101449275362, 50.0414937759336
206.95652173913044, 19.170124481327804
404.05797101449275, 11.203319502074692
559.4202898550725, 6.721991701244804
573.3333333333333, 20
724.0579710144928, 10.871369294605813
906.086956521739, 13.360995850622416];

xdata = data(:,1)/(0.001*519.6*1000);
ydata = data(:,2)/100;

hold on
scatter(xdata,ydata,'k','LineWidth',2)
legend('BRAF = 3nM','BRAF = 10nM','patient data','FontSize',16)
 

ax1=gca;
ax2 = axes('Position', get(ax1, 'Position'),'Color', 'none');
set(ax2, 'XAxisLocation', 'top','YAxisLocation','Right');
set(ax2, 'XLim', get(ax1, 'XLim'),'YLim', get(ax1, 'YLim'));
set(ax2, 'XTick', get(ax1, 'XTick'), 'YTick', 0:100:100);
%Top labels:
%[0:200:1000]/519.6=[0 0.3849 0.7698 1.1547 1.5396 1.9246]
muMdoseLabels = {'\color{blue}0', '\color{blue}0.3849', '\color{blue}0.7698', '\color{blue}1.1547', '\color{blue}1.5396', '\color{blue}1.9246'};
%Dummy labels for right-axis
EmptyTickLabels = {'' '' '' '' ''};
set(ax2, 'XTickLabel', muMdoseLabels ,'YTickLabel',EmptyTickLabels);

xlabel('\color{blue} Dabrafenib concentration (\muM)','FontSize',12)
