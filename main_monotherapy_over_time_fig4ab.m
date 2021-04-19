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

% Set inital conditions
MEK_tot=1.2;
ERK_tot=1.2;
phosph1_tot=0.0003;
phosph2_tot=0.12;
y0=zeros(size(M,1),1);
y0(4)=MEK_tot;
y0(20)=ERK_tot; 
y0(12)=phosph1_tot;
y0(27)=phosph2_tot;


tspan = [0 48*3600];
DBFlist=[0 1 5 10 25 50 100]/10;
TMTlist=DBFlist;%/10;
options = odeset('Mass',M,'MassSingular','yes', 'RelTol',1e-4,'AbsTol',1e-4);
set(gcf,'position',[100,100,1000,600])
newcolors = {'#008000',' #FFC300 ', '#FF5733 ','#C70039','#900C3F','#581845','#000000','#DAF7A6'};

%%%%%%%%%%%%%%%%%%%%%%% Generate FIGURE 4a (dabrafenib monotherapy - ppERK over time).
tiles = tiledlayout(2,2,'TileSpacing','none','Padding','none');
    
BRAF_tot=[0.003,0.003,0.01,0.01];
ATP_tot=[1000, 5000, 1000, 5000]; 
TMT_tot=0;
y0(30)=TMT_tot;
TMT_in=TMT_tot;
    
for k = 1:4
    y0(1)=BRAF_tot(k);
    y0(2)=ATP_tot(k);

    nexttile
    colororder(newcolors)
    for j = 1:length(DBFlist)
        DBF_in=DBFlist(j);
        y0(15)=DBF_in;             
        [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), DBF_in, TMT_in), tspan, y0, options);
        yperk = y(:,22)/1.2;
        ypperk = y(:,26)/1.2;
        hold on
        plot(t,ypperk,'LineWidth',2)
        clear t y
    end
    
    if(k==5)
    legend('0\muM','1\muM', '5\muM', '10\muM', '25\muM','50\muM','100\muM',...
        'Position',[0 100 0 50],'Orientation','horizontal','FontSize',16)
    end

    yticks([0 0.5 1])
    xticks(3600*[0 8 16 24 48])
    xlim([0 tspan(end)])
    xticklabels({'0','8','16','24','48'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)

    xlabel('time (hours)','FontSize',18)
    ylabel('Activated ERK','FontSize',18)
    grid on

end

%%%%%%%%%%%%%%%%%%%%%%% Generate FIGURE 4a (trametinib monotherapy - ppERK over time)

set(gcf,'position',[100,100,1000,600])
newcolors = {'#008000',' #FFC300 ', '#FF5733 ','#C70039','#900C3F','#581845','#000000','#DAF7A6'};

figure
tiles = tiledlayout(2,2,'TileSpacing','none','Padding','none');

BRAF_tot=[0.003,0.003,0.01,0.01];
ATP_tot=[1000, 5000, 1000, 5000]; 
DBF_tot=0;
y0(15)=DBF_tot;
DBF_in=DBF_tot;
    
for k = 1:4
    y0(1)=BRAF_tot(k);
    y0(2)=ATP_tot(k);
      
    nexttile
    colororder(newcolors)
    for j = 1:length(DBFlist)
        TMT_in=TMTlist(j);
        y0(30)=TMT_in;             
        [t,y] = ode15s(@(t,y) mapk_cascade_DAE(y, BRAF_tot(k), ATP_tot(k), DBF_in, TMT_in), tspan, y0, options);
        yperk = y(:,22)/1.2;
        ypperk = y(:,26)/1.2;
        hold on
        plot(t,ypperk,'LineWidth',2)
        
    end

    if(k==5)
    legend('0\muM','0.1\muM', '0.5\muM', '1\muM', '2.5\muM','5\muM','10\muM',...
       'Position',[0 100 0 50],'Orientation','horizontal','FontSize',16)
    end
    yticks([0 0.5 1])
    xticks(3600*[0 8 16 24 48])
    xlim([0 tspan(end)])
    xticklabels({'0','8','16','24','48'})
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'fontsize',14)

    xlabel('time (hours)','FontSize',18)
    ylabel('Activated ERK','FontSize',18)
    grid on

end

