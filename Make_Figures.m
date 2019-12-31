% Script to create figures in the main text and supplement
% This was last updated Dec 30, 2019 by Megan Lickley
% Run this script after running the MainScript_multi.m file for each
% scenario

close all
clear all
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

HomeDir = '/Users/meganlickley/Dropbox (MIT)/Code';
CFC11_filename.LT45 = strcat(HomeDir,'/CFC11/Output/CFC11_fixedLT.mat');
CFC11_filename.LT62p9 = strcat(HomeDir,'/CFC11/Output/CFC11_fixedLT62p9.mat');
CFC11_filename.meanLT = strcat(HomeDir,'/CFC11/Output/CFC11_meanSparcLT.mat');
CFC11_filename.fugitive = strcat(HomeDir,'/CFC11/Output/CFC11_meanSparcLT_fugitive.mat');

CFC12_filename.LT100 = strcat(HomeDir,'/CFC12/Output/CFC12_fixedLT.mat');
CFC12_filename.meanLT = strcat(HomeDir,'/CFC12/Output/CFC12_meanSparcLT.mat');

CFC113_filename.LT85 = strcat(HomeDir,'/CFC113/Output/CFC113_fixedLT.mat');
CFC113_filename.meanLT = strcat(HomeDir,'/CFC113/Output/CFC113_meanSparcLT.mat');
FolderName = 'Figures';


%% Figure 1a: Banks

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

FigHandle = figure(1)
maxylim = 3200;
load(CFC11_filename.meanLT)
p2 = plot([y1:y2],0.001*Banks_TD,'--b','LineWidth',2);
load(CFC11_filename.LT45)
hold on; 
p5 = plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2);
load('Figure_data/WMO2002CFC11data.mat','year','Revised_banks') 
hold on; 
p6 = plot([year(1):year(1)+69],0.001*Revised_banks(1:70),'k','LineWidth',2);
lgd = legend([p2 p5 p6],'SPARC MMM LT','LT = 45yrs','WMO 2003')
lgd.Location = 'northwest';
title('CFC-11 top-down bank estimates'); ylabel('CFC-11 [Gg]'); 
xlim([1955,2016]); ylim([0,maxylim]); xlabel('Year');
set(gca, 'FontSize', 12)
box on

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig1a.pdf');
print(gcf, '-dpdf', str);

%% Figure 1b: Banks and production

clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

FigHandle = figure(2)
maxylim = 3200;
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
load('CFC11/Output/meanvarsCFC11_meanSPARCLT.mat')
Prod_coeff = [1,1.05,1.1];
for ii = 1:3
    Bank(:,ii) = cumsum(Prod_coeff(ii)*Prod_mu-ObsDerivedEmiss);
end
BankSize1 = 5893.9;
Bank = Bank+BankSize1;

pn = plot(years,0.001*Bank', 'LineWidth',2);

load('Figure_data/WMO2002CFC11data.mat','year','Revised_banks') 
hold on; 
p6 = plot([year(1):year(1)+69],0.001*Revised_banks(1:70),'k','LineWidth',2);
set(gca,'FontSize',12); 
lgd = legend([pn; p6],'reported production','+5% reported production','+10% reported production','WMO 2003')
lgd.Location = 'northwest';
xlim([1955,2016]); ylim([0,maxylim]); ylabel('CFC-11 Bank [Gg]'); xlabel('Year')
title('CFC11 top-down bank estimates'); 

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig1b.pdf');
print(gcf, '-dpdf', str);

%% Figure 1c: Bank comparison  
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

FigHandle = figure(3)
maxylim = 3200;
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

load(CFC11_filename.fugitive)
Bank_BPE = Bank(ResampleIndex,:);
y1f = find(years == 1980);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
hold on; 
p3 = boundedline([1980:y2],0.001*MED(y1f:end)',0.001*[LB(y1f:end),UB(y1f:end)],'alpha','cmap',[0.3,0.3,0.3]);
fugemiss = MED;

load(CFC11_filename.meanLT)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,1]);

load(CFC11_filename.LT45)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
hold on; 
p4 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);

load('Figure_data/WMO2002CFC11data.mat','year','Revised_banks') 
hold on; 
p6 = plot([year(1):year(1)+69],0.001*Revised_banks(1:70),'k','LineWidth',2);

% Add TEAP 2009 bank estimates
plot(2008,1420,'+k','MarkerSize',8,'MarkerEdgeColor',[1 0.3 0.3]); %Table 5A-2 from WMO 2010 report
txt = 'TEAP(2009)';
text(2000.5,1358,txt,'Color',[1 0.3 0.3],'FontSize',12)

% Add Ashford estimate
plot(2002,1.68*10^3,'+','MarkerSize',8,'MarkerEdgeColor',[0 0.6 0]); %Page 166 of 
txt = 'Ashford(2004)';
text(1993.5,1.60*10^3,txt,'Color',[0 0.6 0],'FontSize',12)

% Add WMO 2018 estimate
plot([2008,2015],[1420,982],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',10,'MarkerEdgeColor',[0 0 0]);
txt = 'WMO(2018)';
text(2009,900,txt,'Color',[0 0 0],'FontSize',12)

% Add IPCC/TEAP 2019 estimate
plot([2012],[1378],'o','LineWidth',2,'Color',[.8 .2 .8]); 
txt = 'IPCC/TEAP(2019)';
text(2000,3*10^3,txt,'Color',[.8 .2 .8],'FontSize',12)
x = [0.85 0.85];
y = [0.84 0.52];
annotation('textarrow',x,y,'Color',[.8 .2 .8])

lgd = legend([p1 p3 p4 p6],'SPARC MMM LT','SPARC MMM LT, unexpected','LT = 45yrs','WMO 2003')
lgd.Location = 'northwest';

title('CFC-11 BPE and published bank estimates'); ylabel('CFC-11 [Gg]'); 
xlim([1955,2016]); ylim([0,maxylim]); xlabel('Year');
set(gca, 'FontSize', 12)
box on

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig1c.pdf');
print(gcf, '-dpdf', str);

%% Figure 2d: Comparing Banks with constant vs changing lifetime
FigHandle = figure(4)

load(CFC11_filename.meanLT)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,1]);

maxylim = 3200;
hold on; 
p2 = plot([y1:y2],0.001*Banks_TD,'--b','LineWidth',2);

load(CFC11_filename.LT62p9)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
p3 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);
hold on; 
p4 = plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2);

title('CFC-11 BPE bank estimates'); ylabel('CFC-11 [Gg]'); 
xlim([1955,2016]); ylim([0,maxylim]); xlabel('Year')
set(gca, 'FontSize', 12)
box on

lgd = legend([p1 p2 p3 p4],'BPE, SPARC MMM LT','Top-Down SPARC MMM LT','BPE, LT = 62.9yrs','Top-Down, LT = 62.9yrs')
lgd.Location = 'northwest';

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig1d.pdf');
print(gcf, '-dpdf', str);

%% Figure 1e: CFC-12 Banks

FigHandle = figure(5)
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

load(CFC12_filename.meanLT)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,1]);

maxylim = 0.001*1.3*max(mean(Bank_BPE,1));

p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,1]);
hold on; 
p2 = plot([y1:y2],0.001*Banks_TD,'--b','LineWidth',2);


load(CFC12_filename.LT100)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
p3 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);
hold on; 
p4 = plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2);

load('Figure_data/WMO2002CFC12data.mat') 
hold on; 
p5 = plot([year(1):year(end)],0.001*Revised_banks(1:53),'k','LineWidth',2);

% Add WMO 2018 estimate
plot([2008,2015],[394,47],':o','LineWidth',2,'Color',[0 0 0],'MarkerSize',4,'MarkerEdgeColor',[0 0 0]);
txt = 'WMO(2018)';
text(2000,150,txt,'Color',[0 0 0],'FontSize',12)

% ADD TEAP 2009 estimate
plot(2008,394,'+k','MarkerSize',10,'MarkerEdgeColor',[1 0.3 0.3]); %Table 5A-2 from WMO 2010 report
txt = 'TEAP(2009)';
text(2003.5,500,txt,'Color',[1 0.3 0.3],'FontSize',12)

% Add Ashford
plot(2003,0.65*10^3,'+','MarkerSize',8,'MarkerEdgeColor',[0 0.6 0]); %Page 166 of SROC
txt = 'Ashford(2004)';
text(1997.5,0.75*10^3,txt,'Color',[0 0.6 0],'FontSize',12)

title('CFC-12 BPE and published bank estimates'); ylabel('CFC-12 [Gg]');xlabel('Year');
xlim([1955,2016]); ylim([0,1.05*maxylim])
set(gca, 'FontSize', 12)
box on

lgd = legend([p1 p2 p3 p4 p5],'BPE, SPARC MMM LT','Top-Down SPARC MMM LT','BPE, LT = 100yrs','Top-Down, LT = 100yrs','WMO 2003')
lgd.Location = 'northwest';

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig2e.pdf');
print(gcf, '-dpdf', str);

%% Figure 1f: CFC-113 Banks
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

FigHandle = figure(6)
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

load(CFC113_filename.meanLT)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,1]);

maxylim = 0.001*1.6*max(mean(Bank_BPE,1));
hold on; 
p2 = plot([y1:y2],0.001*Banks_TD,'--b','LineWidth',2);

load(CFC113_filename.LT85)
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50);
UB = (prctile(Bank_BPE,97.5)-MED)';
LB = (MED-prctile(Bank_BPE,2.5))';
p3 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);
hold on; 
p4 = plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2);

load('Figure_data/WMO2002CFC113data.mat') 
hold on; 
p5 = plot([year(1):year(end)],0.001*Revised_banks,'k','LineWidth',2);

title('CFC-113 BPE and published bank estimates'); ylabel('CFC-113 [Gg]'); 
xlim([1955,2016]); ylim([0,1.1*maxylim]); xlabel('Year')
set(gca, 'FontSize', 12)
box on

lgd = legend([p1 p2 p3 p4 p5],'BPE, SPARC MMM LT','Top-Down SPARC MMM LT','BPE, LT = 85yrs','Top-Down, LT = 85yrs','WMO 2003');
lgd.Location = 'northwest';

figure_width = 7; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig2f.pdf');
print(gcf, '-dpdf', str);

%% Figure 2: Emissions, production and modeled production
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

FigHandle = figure(7)
clear str Bank_emiss Direct_emiss Total_emiss Prod_prior
str1{1} = CFC11_filename.meanLT ;
str1{2} = CFC12_filename.meanLT ;
str1{3} = CFC113_filename.meanLT ;

for mol_ii = 1:3
    load(str1{mol_ii})
    MED_bank = prctile(Bank(ResampleIndex,:),50);
    MED_emiss = prctile(Emiss(ResampleIndex,:),50);
    MED_prod = prctile(Prod_prior(ResampleIndex,:),50);
    MED_RF = prctile(RF_prior(:,ResampleIndex)',50);
    MED_DE = prctile(DE_prior(:,ResampleIndex)',50);
    Bank_emiss(mol_ii,:) = MED_RF(2:end).*MED_bank(1:end-1);
    Direct_emiss(mol_ii,:) = MED_DE.*MED_prod;
    Total_emiss(mol_ii,:) = MED_emiss;
    Prodvals(mol_ii,:) = MED_prod;
end
CFC11 = [Bank_emiss(1,:)',Direct_emiss(1,2:end)'];
CFC12 = [Bank_emiss(2,:)',Direct_emiss(2,2:end)'];
CFC113 = [Bank_emiss(3,:)',Direct_emiss(3,2:end)'];

subplot(3,1,1)
h1 = area([y1+1:y2],0.001*CFC11);
set(h1(1),'FaceColor',[0.3 0.3 0.3]);
set(h1(2),'FaceColor',[0.6 0.6 0.6]);
FixedLT_line1 = Prodvals(1,:);
hold on; plot([y1:y2],0.001*FixedLT_line1 ,'--r','LineWidth',2)

ylabel('[Gg/yr]');
xlabel('Year'); xlim([y1+1,y2])
set(gca,'FontSize',12); title('SPARC MMM')
legend('Bank Emissions','Direct Emissions','Estimated Production');
title('CFC-11')

subplot(3,1,2)
h1 = area([y1+1:y2],0.001*CFC12);
set(h1(1),'FaceColor',[0.3 0.3 0.3]);
set(h1(2),'FaceColor',[0.6 0.6 0.6]);
FixedLT_line1 = Prodvals(2,:);
hold on; plot([y1:y2],0.001*FixedLT_line1 ,'--r','LineWidth',2)

ylabel('[Gg/yr]');
xlabel('Year'); xlim([y1+1,y2])
set(gca,'FontSize',12); title('SPARC MMM')
title('CFC-12')

subplot(3,1,3)
h1 = area([y1+1:y2],0.001*CFC113);
set(h1(1),'FaceColor',[0.3 0.3 0.3]);
set(h1(2),'FaceColor',[0.6 0.6 0.6]);
FixedLT_line1 = Prodvals(3,:);
hold on; plot([y1:y2],0.001*FixedLT_line1 ,'--r','LineWidth',2)

ylabel('[Gg/yr]');
xlabel('Year'); xlim([y1+1,y2])
set(gca,'FontSize',12); title('SPARC MMM')
title('CFC-113')

figure_width = 8; % in inches
figure_height = 12; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig2.pdf');
print(gcf, '-dpdf', str);
%% Figure 3:  Emissions posterior and obs-derived emissions

clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

SCEN_FILE{1} = CFC11_filename.meanLT;
SCEN_FILE{2} = CFC12_filename.meanLT;
SCEN_FILE{3} = CFC113_filename.meanLT;

Mole_Name{1} = 'CFC-11';
Mole_Name{2} = 'CFC-12';
Mole_Name{3} = 'CFC-113';

Input_FileName{1} = 'CFC11';
Input_FileName{2} = 'CFC12';
Input_FileName{3} = 'CFC113';

LT_const(1) = 62.9;
LT_const(2) = 100;
LT_const(3) = 85;

ppt_to_tonnes(1) =22602.38457; % Conversion of ppt to tonnes of CFC-11
ppt_to_tonnes(2) =19895.36; % Conversion of ppt to tonnes of CFC-12
ppt_to_tonnes(3) =30834.68811; % Conversion of ppt to tonnes of CFC-113

FigHandle = figure(8);
for scen_ii = 1:3
    load(SCEN_FILE{scen_ii})
    subplot(1,3,scen_ii)
    MED = prctile(Emiss(ResampleIndex,:),50);
    UB = (prctile(Emiss(ResampleIndex,:),97.5)-MED);
    LB = (MED-prctile(Emiss(ResampleIndex,:),2.5));
    p1 = boundedline([y1:y2],0.001*MED',[0.001*LB',0.001*UB'],'alpha','cmap',[0.3,0.3,0.3]);
    UBmax = (prctile(Emiss(ResampleIndex,:),99.5));
    LBmin = (prctile(Emiss(ResampleIndex,:),0.5));
    hold on; plot([y1:y2],0.001*UBmax,'--k');
    hold on; plot([y1:y2],0.001*LBmin,'--k');
    hold on; 
    p2 = plot([y1:y2],0.001*ObsDerivedEmiss,'r','LineWidth',2);
    str = strcat(Mole_Name{scen_ii}, ' Emissions');
    title(str);
    str = strcat(Mole_Name{scen_ii}, ' emissions [Gg yr^{-1}]');
    box on; ylabel(str); xlabel('year'); ylim([0,0.001*1.7*max(ObsDerivedEmiss)])
    set(gca,'FontSize',12)
    xlim([1955,2016])
    
    str = strcat(Input_FileName{scen_ii},'/Input/wmo2018.mat');
    load(str);
    ytmp1 = find(wmo_yr == y1);
    ytmp2 = find(wmo_yr == y2);
    conc_obs = wmo_conc(ytmp1:ytmp2+1);
    Emiss_constantLT = ppt_to_tonnes(scen_ii)*(conc_obs(2:end)'-conc_obs(1:end-1)'.*exp(-1./LT_const(scen_ii)'));
    hold on; 
    p3 = plot([y1:y2],0.001*Emiss_constantLT,'g','LineWidth',2);
    str = strcat('LT =',' ',num2str(LT_const(scen_ii)),' yrs');
    legend([p1 p2 p3],{'bayesian posterior','LT = MMM', str})
end

figure_width = 20; % in inches
figure_height = 5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig3.pdf');
print(gcf, '-dpdf', str);


%% Fig 3 insert: Prior Emissions and obs-derived emissions insert

FigHandle = figure(9);
for scen_ii = 1:3
    load(SCEN_FILE{scen_ii})
    subplot(1,3,scen_ii)
    MED = prctile(Emiss(ResampleIndex,:),50);
    UB = (prctile(Emiss(ResampleIndex,:),97.5)-MED);
    LB = (MED-prctile(Emiss(ResampleIndex,:),2.5));
    p1 = boundedline([y1:y2],0.001*MED',[0.001*LB',0.001*UB'],'alpha','cmap',[0.3,0.3,0.3]);
    UBmax = (prctile(Emiss(ResampleIndex,:),99.5));
    LBmin = (prctile(Emiss(ResampleIndex,:),0.5));
    hold on; plot([y1:y2],0.001*UBmax,'--k');
    hold on; plot([y1:y2],0.001*LBmin,'--k');
    hold on; 
    p2 = plot([y1:y2],0.001*ObsDerivedEmiss,'r','LineWidth',2);
    
    str = strcat(Mole_Name{scen_ii}, ' emissions [Gg yr^{-1}]');
    box on; ylabel(str); xlabel('year'); ylim([0,0.001*0.3*max(ObsDerivedEmiss)])
    set(gca,'FontSize',12)
    xlim([2010,2016])
    
    str = strcat(Input_FileName{scen_ii},'/Input/wmo2018.mat');
    load(str);
    ytmp1 = find(wmo_yr == y1);
    ytmp2 = find(wmo_yr == y2);
    conc_obs = wmo_conc(ytmp1:ytmp2+1);
    Emiss_constantLT = ppt_to_tonnes(scen_ii)*(conc_obs(2:end)'-conc_obs(1:end-1)'.*exp(-1./LT_const(scen_ii)'));
    hold on; 
    p3 = plot([y1:y2],0.001*Emiss_constantLT,'g','LineWidth',2);
end

figure_width = 10; % in inches
figure_height = 2.5; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);
str = strcat(FolderName,'/Fig3insert.pdf');
print(gcf, '-dpdf', str);

 
%% Figure 4 EESC curve
% The equation for each species is as follows:
% (no. atoms)*MMR*(air/species)*(fractional halogen release)*(Br scaling factor if needed)
% Fractional halon release and Br scaling factor obtained from: 
% Newman et al., (2007) https://doi.org/10.5194/acp-7-4537-2007
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

FigHandle = figure(10)
fugitive = 0;
if fugitive
    str1{1} = CFC11_filename.fugitive;
else
    str1{1} = CFC11_filename.meanLT;
end
str1{2} = CFC12_filename.meanLT ;
str1{3} = CFC113_filename.meanLT ;

if fugitive 
    str2{1} = 'CFC11/Output/meanvarsCFC11_meanSparcLT_fugitive.mat'; 
else
    str2{1} = 'CFC11/Output/meanvarsCFC11_meanSparcLT.mat';
end
str2{2} = 'CFC12/Output/meanvarsCFC12_meanSparcLT.mat';
str2{3} = 'CFC113/Output/meanvarsCFC113_meanSparcLT.mat';

y1 = 1955; % beginning of BPE simulation period
y2 = 2016; % end of BPE simulation period
y3 = 2100; % end year of forward simulation
y4 = 2020; %Scen 2, destroy bank in 2020
y5 = 2000; %Scen 3, opportunity lost for not destroying banks in 2000
years = y1:y3;

ni = [3,2,3]; % number of chlorine atoms
air = 28.97; %(g mol-1)
fi = [0.99, 0.86, 0.90]; % From Newman et al., 2007 A new formulation of equivalent effective stratospheric chlorine
mole_weight = [137.37,120.91,187.37];

ppt_to_tonnes(1) = 22602.38457; % Conversion for CFC11
ppt_to_tonnes(2) = 19895.36;  % Conversion for CFC12
ppt_to_tonnes(3) = 30834.68811; % Conversion for CFC113

% start at last time period of Bayesian analysis and run the model forward
% for three scenarios: 
%1. Banks are not destroyed (Business as usual)
%2. Banks are destroyed in 2020
%3. Banks were destroyed in 2000
load('Figure_data/WMO2018data.mat')
yindx1 = find(yr_wmo == y1);
yindx2 = find(yr_wmo == y2);
yindx4 = find(yr_wmo == y4);

for mol_ii = 1:3
    
    load(str1{mol_ii})
    load(str2{mol_ii})
    MMM_Bank = prctile(Bank(ResampleIndex,:),50);
    MMM_RF = prctile(RF_prior(:,ResampleIndex)',50);
    MMM_LT = LT;
    MMM_Emiss = prctile(Emiss(ResampleIndex,:),50);
    clear Bank
    
    ytmp1 = find(years == y1);
    ytmp2 = find(years == y2);
    ytmp3 = find(years == y3);
    ytmp4 = find(years == y4);
    ytmp5 = find(years == y5);
    %Scenario 1. Banks are not destroyed
   
    % Concentrations as in WMO2018 up until 2016
    Conc_scen1(ytmp1:ytmp2,mol_ii) = WM02018(yindx1:yindx2,mol_ii); 
    
    % Estimate the concentration going forward with Bank release
    Bank(ytmp2) = MMM_Bank(end);
    RF = MMM_RF(end); % Keep RF constant going forward
    emissions(ytmp2) = MMM_Emiss(end);
    Emission_timeseries(ytmp2-56:ytmp2,mol_ii) = MMM_Emiss(end-56:end);
    for t = ytmp2+1:ytmp3
        Emission_timeseries(t,mol_ii) = RF*Bank(t-1);
        Bank(t) = (1-RF)*Bank(t-1);
        Conc_scen1(t,mol_ii) = Conc_scen1(t-1,mol_ii)*exp(-1/MMM_LT(end))+Emission_timeseries(t-1,mol_ii)/ppt_to_tonnes(mol_ii);
    end
    
    %Scenario 2. Banks are destroyed in 2020
    
    % Concentrations are the same as Scen 1 until 2020 (ytmp4)
    
    Conc_scen2(ytmp1:ytmp4,mol_ii) = Conc_scen1(ytmp1:ytmp4,mol_ii);
    for t = ytmp4+1:ytmp3
        Conc_scen2(t,mol_ii) = Conc_scen2(t-1,mol_ii)*exp(-1/MMM_LT(end));
    end
    
    %Scenario 3. Banks were destroyed in 2000 and no further production
    % Concentrations are the same as scen 1 until the year 2000
    Conc_scen3(ytmp1:ytmp5,mol_ii) = Conc_scen1(ytmp1:ytmp5,mol_ii);
    for t = ytmp5+1:ytmp3
        Conc_scen3(t,mol_ii) = Conc_scen3(t-1,mol_ii)*exp(-1/MMM_LT(end-(ytmp2-ytmp5)+1));
    end    
end

% For Table 1 values: 
Scen1_emiss = sum(Emission_timeseries(ytmp5:ytmp3,:),1);
Scen2_emiss = sum(Emission_timeseries(ytmp5:ytmp4,:),1);

Scen1_emiss
Scen2_emiss

GWP = [4660,10200,5820];
Scen1_co2eq = Scen1_emiss.*GWP;
Scen2_co2eq = Scen2_emiss.*GWP;

Scen1_co2eq
Scen2_co2eq

% Estimating EESC using age of air spectrum of gamma = 5.5 yrs
alpha = 65; %For polar latitudes
en = 0;
air = 28.97; %(g mol-1)

% Mean release times
%Mean release times for mean age of air of 5.5 years from Table 2 of Engel et al. 2018 (https://www.atmos-chem-phys.net/18/601/2018/acp-18-601-2018.pdf)
% WMO 2018 uses mrtd values as follows, here we use 5.5 following Newman
% 2007. 

mrtd(1) = 5.5;     %CFC-11 
mrtd(2) = 5.9;     %CFC-12                
mrtd(3) = 5.8;     %CFC-113             
mrtd(4) = 8.3;     %CFC-114               
mrtd(5) = 10.1;     %CFC-115              
mrtd(6) = 5.5;    %Carbon tetrachlorine
mrtd(7) = 5.6;     %Methyl chloroform
mrtd(8) = 7.0;     %HCFC-22 
mrtd(9) = 5.8;     %HCFC-141b
mrtd(10) = 6.5;     %HCFC-142b
mrtd(11) = 5.5;     %Halon-1211
mrtd(12) = 5.5;     %Halon-1202
mrtd(13) = 6.2;     %Halon-1301
mrtd(14) = 5.5;     %Halon-2402
mrtd(15) = 5.5;     %Methyl bromide
mrtd(16) = 5.8;     %Methyl chloride

ODS_data = WM02018(:,1:16);

% Construct inverse Gaussian function
% constructing the age spectrum as an inverse gaussian function following Eq. 9 from Waugh and Hall, (2002) (https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000RG000101)
% instead of age of air, we use the mean relese times for each molecule
% separately.

% Here, I am using Engel et al.'s delta caculation of lambda = delta^2/mrtd = .7
% this implies delta = sqrt(aoa*.7)

% Newman et al. uses the delta definition of half the aoa. This is also
% used by WMO 2018. In that case delta = 5.5/2 for Newman and mrtd/2 for
% WMO (2018) (although they were not explicit, they may have also used
% 5.5/2)

inversegaussianfunction = zeros(length(yr_wmo),size(mrtd,2));
for j = 1:size(mrtd,2)
    
    %aoa = mrtd(j);
    aoa = 5.5;
    delta = aoa/2;
    
    count = 1;
    for t = 0:length(yr_wmo)-1
        inversegaussianfunction(count,j) = 1/(2*delta*sqrt(pi*(t/aoa)^3))*exp(-aoa^2*(t/aoa-1)^2./(4*delta^2*(t/aoa)));
        count = count+1;
    end
end

% fractional release values: (Newman et al. (2007))
fr(1) = .99;    %CFC-11                 .99 (Newman)
fr(2) = .86;    %CFC-12                 .86 (Newman)
fr(3) = .90;    %CFC-113                .90 (Newman)
fr(4) = .40;    %CFC-114                .40 (Newman)
fr(5) = .15;     %CFC-115                .15 (Newman)
fr(6) = 1;      %Carbon tetrachlorine   1   (Newman)
fr(7) = .99;    %Methyl chloroform      .99 (Newman)
fr(8) = .41;    %HCFC-22                .41 (Newman)
fr(9) = .90;    %HCFC-141b              .90 (Newman)
fr(10) = .65;   %HCFC-142b              .29 (Newman)
fr(11) = 1;     %Halon-1211             1   (Newman)
fr(12) = 1;     %Halon-1202             1   (Newman)
fr(13) = .8;   %Halon-1301             .8  (Newman)
fr(14) = 1;     %Halon-2402             1   (Newman)
fr(15) = .99;   %Methyl bromide         .99 (Newman)
fr(16) = .91;   %Methyl chloride        .91 (Newman)

clearvars ODS_data_EESC;
ODS_data_EESC = zeros(length(yr_wmo),size(mrtd,2));

for i = 1:length(yr_wmo)
    if i <= 20
        endint = 1;
        endint2 = 1:i;
    else
        endint = i-19;
        endint2 = 1:20;
    end
    ODS_data_EESC(i,1) = 3*fr(1)*nansum(ODS_data(endint:i,1).*flipud(inversegaussianfunction(endint2,1)));              %CFC-11                 (CCl3F) *3 (air./137.37)
    ODS_data_EESC(i,2) = 2*fr(2)*nansum(ODS_data(endint:i,2).*flipud(inversegaussianfunction(endint2,2)));              %CFC-12                 (CCl2F2) *2 (air./120.91)
    ODS_data_EESC(i,3) = 3*fr(3)*nansum(ODS_data(endint:i,3).*flipud(inversegaussianfunction(endint2,3)));              %CFC-113                (CCl2FCClF2) *3 (air./187.37)
    ODS_data_EESC(i,4) = 2*fr(4)*nansum(ODS_data(endint:i,4).*flipud(inversegaussianfunction(endint2,4)));              %CFC-114                (CClF2CClF2) (air./170.92)
    ODS_data_EESC(i,5) = 1*fr(5)*nansum(ODS_data(endint:i,5).*flipud(inversegaussianfunction(endint2,5)));              %CFC-115                (CClF2CF3) (air./154.47)
    ODS_data_EESC(i,6) = 4*fr(6)*nansum(ODS_data(endint:i,6).*flipud(inversegaussianfunction(endint2,6)));              %Carbon tetrachlorine   (CCl4) (air./153.82)
    ODS_data_EESC(i,7) = 3*fr(7)*nansum(ODS_data(endint:i,7).*flipud(inversegaussianfunction(endint2,7)));              %Methyl chloroform      (CH3CCl3) (air./133.4)
    ODS_data_EESC(i,8) = 1*fr(8)*nansum(ODS_data(endint:i,8).*flipud(inversegaussianfunction(endint2,8)));              %HCFC-22                (CHClF3) (air./86.47)
    ODS_data_EESC(i,9) = 2*fr(9)*nansum(ODS_data(endint:i,9).*flipud(inversegaussianfunction(endint2,9)));              %HCFC-141b              (CH3CCl2F) (air./116.95)
    ODS_data_EESC(i,10) = 1*fr(10)*nansum(ODS_data(endint:i,10).*flipud(inversegaussianfunction(endint2,10)));          %HCFC-142b              (CH3CClF2) (air./100.5)
    ODS_data_EESC(i,11) = 1*fr(11)*nansum(ODS_data(endint:i,11).*flipud(inversegaussianfunction(endint2,11)));          %Halon-1211             (CBrClF2) %for Cl (air./165.36)
    
    ODS_data_EESC(i,12) = 1*fr(11)*alpha*nansum(ODS_data(endint:i,11).*flipud(inversegaussianfunction(endint2,11)));    %Halon-1211             (CBrClF2) %for Br (air./165.36)
    ODS_data_EESC(i,13) = 2*fr(12)*alpha*nansum(ODS_data(endint:i,12).*flipud(inversegaussianfunction(endint2,12)));    %Halon-1202             (CBr2F2) (air./209.82)
    ODS_data_EESC(i,14) = 1*fr(13)*alpha*nansum(ODS_data(endint:i,13).*flipud(inversegaussianfunction(endint2,13)));    %Halon-1301             (CBrF3) (air./148.91)
    ODS_data_EESC(i,15) = 2*fr(14)*alpha*nansum(ODS_data(endint:i,14).*flipud(inversegaussianfunction(endint2,14)));    %Halon-2402             (CBrF2CBrF2) (air./259.8)
    ODS_data_EESC(i,16) = 1*fr(15)*alpha*nansum(ODS_data(endint:i,15).*flipud(inversegaussianfunction(endint2,15)));    %Methyl bromide         (CH3Br) (air./94.94)
   
    ODS_data_EESC(i,17) = 1*fr(16)*nansum(ODS_data(endint:i,16).*flipud(inversegaussianfunction(endint2,16)));          %Methyl chloride        (CH3Cl) (air./50.49)
    

    EESC_scen1(i,1) = 3*fr(1)*nansum(Conc_scen1(endint:i,1).*flipud(inversegaussianfunction(endint2,1)));              %CFC-11                 (CCl3F) *3 (air./137.37)
    EESC_scen1(i,2) = 2*fr(2)*nansum(Conc_scen1(endint:i,2).*flipud(inversegaussianfunction(endint2,2)));              %CFC-12                 (CCl2F2) *2 (air./120.91)
    EESC_scen1(i,3) = 3*fr(3)*nansum(Conc_scen1(endint:i,3).*flipud(inversegaussianfunction(endint2,3)));              %CFC-113                (CCl2FCClF2) *3 (air./187.37)
    
    EESC_scen2(i,1) = 3*fr(1)*nansum(Conc_scen2(endint:i,1).*flipud(inversegaussianfunction(endint2,1)));              %CFC-11                 (CCl3F) *3 (air./137.37)
    EESC_scen2(i,2) = 2*fr(2)*nansum(Conc_scen2(endint:i,2).*flipud(inversegaussianfunction(endint2,2)));              %CFC-12                 (CCl2F2) *2 (air./120.91)
    EESC_scen2(i,3) = 3*fr(3)*nansum(Conc_scen2(endint:i,3).*flipud(inversegaussianfunction(endint2,3)));              %CFC-113                (CCl2FCClF2) *3 (air./187.37)

    EESC_scen3(i,1) = 3*fr(1)*nansum(Conc_scen3(endint:i,1).*flipud(inversegaussianfunction(endint2,1)));              %CFC-11                 (CCl3F) *3 (air./137.37)
    EESC_scen3(i,2) = 2*fr(2)*nansum(Conc_scen3(endint:i,2).*flipud(inversegaussianfunction(endint2,2)));              %CFC-12                 (CCl2F2) *2 (air./120.91)
    EESC_scen3(i,3) = 3*fr(3)*nansum(Conc_scen3(endint:i,3).*flipud(inversegaussianfunction(endint2,3)));              %CFC-113                (CCl2FCClF2) *3 (air./187.37)

end


EESC_rest = sum(ODS_data_EESC(:,4:end),2);
EESC_tot1 = 0.001*EESC_rest+0.001*sum(EESC_scen1,2);
EESC_tot2 = 0.001*EESC_rest+0.001*sum(EESC_scen2,2);
EESC_tot3 = 0.001*EESC_rest+0.001*sum(EESC_scen3,2);
EESC_WMO = 0.001*sum(ODS_data_EESC(:,1:end),2);

Panel1_labels{1} = 'CH3Cl';
Panel1_labels{2} = 'CCl4';
Panel1_labels{3} = 'CFCs';
Panel1_labels{4} = 'HCFCs';
Panel1_labels{5} = 'Halons';
Panel1_labels{6} = 'CH3Br';
Panel1_labels{7} = 'CH3CCl3';

EESC_panel1(:,1) = 0.001*ODS_data_EESC(:,17);
EESC_panel1(:,2) = 0.001*ODS_data_EESC(:,6);
EESC_panel1(:,3) = 0.001*sum(ODS_data_EESC(:,1:5),2);
EESC_panel1(:,4) = 0.001*sum(ODS_data_EESC(:,8:10),2);
EESC_panel1(:,5) = 0.001*sum(ODS_data_EESC(:,11:15),2);
EESC_panel1(:,6) = 0.001*ODS_data_EESC(:,16);
EESC_panel1(:,7) = 0.001*ODS_data_EESC(:,7);

Panel2_labels{1} = 'CH3Cl';
Panel2_labels{2} = 'CCl4';
Panel2_labels{3} = 'CH3Br';
Panel2_labels{4} = 'CH3CCl3';
Panel2_labels{5} = 'Halons';
Panel2_labels{6} = 'HCFCs';
Panel2_labels{7} = 'CFCs';
Panel2_labels{8} = 'Scenario 1';
Panel2_labels{9} = 'Scenario 2';
Panel2_labels{10} = 'Scenario 3';

EESC_panel2(:,1) = 0.001*ODS_data_EESC(:,17);
EESC_panel2(:,2) = 0.001*ODS_data_EESC(:,6);
EESC_panel2(:,7) = 0.001*sum(ODS_data_EESC(:,1:5),2);
EESC_panel2(:,6) = 0.001*sum(ODS_data_EESC(:,8:10),2);
EESC_panel2(:,5) = 0.001*sum(ODS_data_EESC(:,11:15),2);
EESC_panel2(:,3) = 0.001*ODS_data_EESC(:,16);
EESC_panel2(:,4) = 0.001*ODS_data_EESC(:,7);


subplot(1,2,1)
h = area(yr_wmo',EESC_panel1); xlim([1975,2099]); ylabel('EESC ppbv'); ylim([0,4.5]);
h(1).FaceColor = [0 0 1];
h(2).FaceColor = [0.5 0 1];
h(3).FaceColor = [1 0 1];
h(4).FaceColor = [0.5 0.5 0.5];
h(5).FaceColor = [1 0 0];
h(6).FaceColor = [1 0.5 0];
h(7).FaceColor = [1 1 0];
%legend(Panel1_labels)
set(gca,'FontSize',12); title('EESC')
grid on

subplot(1,2,2)
h = area(yr_wmo',EESC_panel2); xlim([1975,2099]); ylabel('EESC ppbv'); ylim([0,4.5]);
h(1).FaceColor = [0 0 1];
h(2).FaceColor = [0.5 0 1];
h(7).FaceColor = [1 0 1];
h(6).FaceColor = [0.5 0.5 0.5];
h(5).FaceColor = [1 0 0];
h(3).FaceColor = [1 0.5 0];
h(4).FaceColor = [1 1 0];
hold on; plot(yr_wmo',EESC_tot1,':k','LineWidth',2)
hold on; plot(yr_wmo',EESC_tot2,'--k','LineWidth',2)
hold on; plot(yr_wmo',EESC_tot3,'-.r','LineWidth',2)
legend(Panel2_labels)
set(gca,'FontSize',12); title('EESC')
grid on


yind = find(yr_wmo == 1980);
EESC1980 =  EESC_tot1(yind);

xind1 = find(sum(EESC_tot1(36:end,:),2)<EESC1980,1)+35;
xind2 = find(sum(EESC_tot2(36:end,:),2)<EESC1980,1)+35;
xind3 = find(sum(EESC_tot3(36:end,:),2)<EESC1980,1)+35;
xind4 = find(sum(EESC_WMO(36:end,:),2)<EESC1980,1)+35;
yr_wmo(xind4)

hold on; plot([1980,1980],[0,EESC1980],'--k','LineWidth',2)
hold on; plot([1980,yr_wmo(xind1)],[EESC1980,EESC1980],'--k','LineWidth',2)
hold on; plot([yr_wmo(xind1),yr_wmo(xind1)],[0,EESC1980],'--k','LineWidth',2)
hold on; plot([yr_wmo(xind2),yr_wmo(xind2)],[0,EESC1980],'--k','LineWidth',2)
hold on; plot([yr_wmo(xind3),yr_wmo(xind3)],[0,EESC1980],'--k','LineWidth',2)


figure_width = 12; % in inches
figure_height = 6; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig4.pdf');
print(gcf, '-dpdf', str);


%%  Fig 5: concentrations with uncertainties for three scenarios

clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName


y1 = 1960;
y2 = 2016;
y3 = 2100;
y4 = 2019;
y5 = 2029;
years = y1:y3;

ppt_to_tonnes(1) = 22602.38457; %CFC11
ppt_to_tonnes(2) = 19895.36;  % CFC12
ppt_to_tonnes(3) = 30834.68811; %CFC113

% start at last time period of Bayesian analysis and run the model forward
% for three scenarios: 
%1. Concentrations with extra production
%2. Concentrations without extra production
%3. WMO2018
load('Figure_data/WMO2018data.mat')
yindx1 = find(yr_wmo == y1);
yindx2 = find(yr_wmo == y2);
yindx3 = find(yr_wmo == y3);
yindx4 = find(yr_wmo == y4);
yindx5 = find(yr_wmo == y5);

%%%%%%%%%%%%%%%%%%%%%%%
% CFC-11 scenarios:
mol_ii = 1;
scen_ii = 2;
load(CFC11_filename.meanLT)
load('CFC11/Output/meanvarsCFC11_meanSparcLT.mat')

ytmp1 = find(years == y1);
ytmp2 = find(years == y2);
ytmp3 = find(years == y3);
ytmp4 = find(years == y4);
ytmp5 = find(years == y5);

%i) WMO 2018 projection

Conc(ytmp1:ytmp3,1,mol_ii) = WM02018(yindx1:yindx3,mol_ii); 

%ii) * Bayesian projection assuming no fugitive production
Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2016

Bank_samps = Bank(ResampleIndex,end);
RF_samps = RF_prior(end,ResampleIndex)';
Emiss_samps = Emiss(ResampleIndex,end);
DE_samps = DE_prior(end,ResampleIndex)';

clear Bank RF Emiss Prod_prior DE_prior
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);

tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    tind = tind+1;
end
for t = ytmp2+3:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end

Conc(ytmp2+1:ytmp3,scen_ii,mol_ii) = prctile(Conc_samps(:,ytmp2+1:ytmp3),50);
Conc_UB(:,scen_ii,mol_ii) = prctile(Conc_samps,84)-prctile(Conc_samps,50);
Conc_LB(:,scen_ii,mol_ii) = prctile(Conc_samps,50)-prctile(Conc_samps,16);

%iii) Bayesian projection assuming fugitive production that stops in 2019
scen_ii = 3;
load(CFC11_filename.fugitive);
load('CFC11/Output/meanvarsCFC11_meanSparcLT_fugitive.mat'); 
Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2016

Bank_samps = Bank(ResampleIndex,end);
RF_samps = RF_prior(end,ResampleIndex)';
Emiss_samps = Emiss(ResampleIndex,end);
Prod_samps = Prod_prior(ResampleIndex,end);
DE_samps = DE_prior(end,ResampleIndex)';

clear Bank RF Emiss Prod_prior DE_prior
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);

tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    tind = tind+1;
end

for t = ytmp2+3:ytmp4
    emissions(:,t) = RF_samps.*Bank_samps(:,tind)+DE_samps.*Prod_samps;
    Direct_emiss(t,scen_ii,mol_ii) = median(DE_samps.*Prod_samps);
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
for t = ytmp4+1:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Direct_emiss(t,scen_ii,mol_ii) = 0;
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end

Conc(ytmp2+1:ytmp3,scen_ii,mol_ii) = prctile(Conc_samps(:,ytmp2+1:ytmp3),50);
Conc_UB(:,scen_ii,mol_ii) = prctile(Conc_samps,84)-prctile(Conc_samps,50);
Conc_LB(:,scen_ii,mol_ii) = prctile(Conc_samps,50)-prctile(Conc_samps,16);


%iv) Bayesian projection assuming fugitive production stops in 2029  (If
%this scenario then I would use the bayesian estimate of what the fugitive
%production was and also run the bank scenario another 10 years to be
%consistent)
scen_ii = 4;
load(CFC11_filename.fugitive);
load('CFC11/Output/meanvarsCFC11_meanSparcLT_fugitive.mat'); 

Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2016

Bank_samps = Bank(ResampleIndex,end);
RF_samps = RF_prior(end,ResampleIndex)';
Emiss_samps = Emiss(ResampleIndex,end);
Prod_samps = Prod_prior(ResampleIndex,end);
DE_samps = DE_prior(end,ResampleIndex)';

clear Bank RF Emiss Prod_prior DE_prior
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);

tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    tind = tind+1;
end

for t = ytmp2+3:ytmp5
    emissions(:,t) = RF_samps.*Bank_samps(:,tind)+DE_samps.*Prod_samps;
    Direct_emiss(t,scen_ii,mol_ii) = median(DE_samps.*Prod_samps);
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
for t = ytmp5+1:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Direct_emiss(t,scen_ii,mol_ii) = 0;
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end

Conc(ytmp2+1:ytmp3,scen_ii,mol_ii) = prctile(Conc_samps(:,ytmp2+1:ytmp3),50);
Conc_UB(:,scen_ii,mol_ii) = prctile(Conc_samps,51)-prctile(Conc_samps,50);
Conc_LB(:,scen_ii,mol_ii) = prctile(Conc_samps,50)-prctile(Conc_samps,49);

FigHandle = figure(12); 
subplot(1,3,1)
scen_ii = 2;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.2,0.2,0.2]);
hold on; 
p1 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'Color',[0.2,0.2,0.2],'LineWidth',2)
scen_ii = 3;
hold on; 
p2 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2);
%p2 = boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0,0,0.7]);
hold on; 
scen_ii = 4;
p3 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2);
hold on; 
p4 = plot([y1:y3],Conc(:,1,mol_ii),'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
 ylim([0,350]);
ylabel('CFC-11 [ppt]')

box on; 
legend([p1 p2 p3 p4], 'Scen 1', 'Scen 2', 'Scen 3', 'WMO 2018')
set(gca,'FontSize',12); 


%%%%%%%%%%%%%%%%%%%%
% CFC-12 scenarios
%i) WMO 2018 projection
%ii) * Bayesian projection

mol_ii = 2;
str1 = CFC12_filename.meanLT;
str2 = 'CFC12/Output/meanvarsCFC12_meanSparcLT.mat';

load(str1)
load(str2)
        
ytmp2 = find(years == y2);
ytmp1 = find(years == y1);
ytmp3 = find(years == y3);

%Scen 1: WMO 2018 projection
Conc(ytmp1:ytmp3,1,mol_ii) = WM02018(yindx1:yindx3,mol_ii); 
% All other scenarios start with WMO2018
scen_ii = 2;
Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2018

% Estimate Emissions from banks and Emissions from non banks for
% all posterior samples.  Take the median and 1-sigma estimates at
% the end

Bank_samps = Bank(ResampleIndex,end);
RF_samps = RF_prior(end,ResampleIndex)';
Emiss_samps = Emiss(ResampleIndex,end);
Prod_samps = Prod_prior(ResampleIndex,end);
DE_samps = DE_prior(end,ResampleIndex)';

clear Bank RF_prior Emiss Prod_prior DE_prior
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);
tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    tind = tind+1;
end
for t = ytmp2+3:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end

Conc(ytmp2+1:ytmp3,scen_ii,mol_ii) = prctile(Conc_samps(:,ytmp2+1:ytmp3),50);
Conc_UB(:,scen_ii,mol_ii) = prctile(Conc_samps,82)-prctile(Conc_samps,50);
Conc_LB(:,scen_ii,mol_ii) = prctile(Conc_samps,50)-prctile(Conc_samps,16);

subplot(1,3,2)
scen_ii = 2;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.2,0.2,0.2]);
hold on;
p1 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2,'Color',[0.2,0.2,0.2]);
hold on;
p2 = plot([y1:y3],Conc(:,1,mol_ii),'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
ylabel('CFC-12 [ppt]'); ylim([0,750]);

box on; 
legend([p1 p2], 'Scen 1', 'WMO 2018')
set(gca,'FontSize',12); 

%%%%%%%%%%%%%%%%%
% CFC-113 scenarios
mol_ii = 3;
str1 = CFC113_filename.meanLT;
str2 = 'CFC113/Output/meanvarsCFC113_meanSparcLT.mat';

%i) WMO 2018 projection
load(str1)
load(str2)
        
ytmp2 = find(years == y2);
ytmp1 = find(years == y1);
ytmp3 = find(years == y3);

Conc(ytmp1:ytmp3,1,mol_ii) = WM02018(yindx1:yindx3,mol_ii); 

%ii) *Bayesian projection assuming no additional emissions past 2018
scen_ii = 2;
Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2016

Bank_samps = Bank(ResampleIndex,end);
RF_samps = RF_prior(end,ResampleIndex)';
Emiss_samps = Emiss(ResampleIndex,end);
Prod_samps = Prod_prior(ResampleIndex,end);
DE_samps = DE_prior(end,ResampleIndex)';

clear Bank RF Emiss Prod_prior DE_prior
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);

tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    tind = tind+1;
end
for t = ytmp2+3:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind)+DE_samps.*Prod_samps;
    Direct_emiss(t,scen_ii,mol_ii) = median(DE_samps.*Prod_samps);
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
Conc(ytmp2+1:ytmp3,scen_ii,mol_ii) = prctile(Conc_samps(:,ytmp2+1:ytmp3),50);
Conc_UB(:,scen_ii,mol_ii) = prctile(Conc_samps,51)-prctile(Conc_samps,50);
Conc_LB(:,scen_ii,mol_ii) = prctile(Conc_samps,50)-prctile(Conc_samps,49);

%iii) **Bayesian projection assuming 7.2 Gg/year out to 2029

scen_ii = 3;
Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2016

Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);

tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    tind = tind+1;
end

for t = ytmp2+3:ytmp5
    emissions(:,t) = RF_samps.*Bank_samps(:,tind)+7.2*10^3;
    Bank_emiss(t,scen_ii,mol_ii) = median(RF_samps.*Bank_samps(:,tind));
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
for t = ytmp5+1:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end

Conc(ytmp2+1:ytmp3,scen_ii,mol_ii) = prctile(Conc_samps(:,ytmp2+1:ytmp3),50);

% Upper Bound (7.2 + 5 Gg)
Conc(ytmp1:ytmp2+2,scen_ii,mol_ii) = WM02018(yindx1:yindx2+2,mol_ii); %2016
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);

tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    tind = tind+1;
end

for t = ytmp2+3:ytmp5
    emissions(:,t) = RF_samps.*Bank_samps(:,tind)+(7.2+5)*10^3;
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
for t = ytmp5+1:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
Conc_UB(:,scen_ii,mol_ii) = prctile(Conc_samps,50)'-Conc(:,scen_ii,mol_ii);

% Lower Bound (7.2 - 5 Gg)
Conc_samps = zeros(size(ResampleIndex,1),ytmp3);
Conc_samps(:,ytmp1:ytmp2+2) = repmat(Conc(ytmp1:ytmp2+2,scen_ii,mol_ii)',size(ResampleIndex,1),1);
tind = 1;
for t = ytmp2+1:ytmp2+2
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind)+(1-DE_samps).*Prod_samps;
    tind = tind+1;
end

for t = ytmp2+3:ytmp5
    emissions(:,t) = RF_samps.*Bank_samps(:,tind)+(7.2-5)*10^3;
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
for t = ytmp5+1:ytmp3
    emissions(:,t) = RF_samps.*Bank_samps(:,tind);
    Bank_samps(:,tind+1) = (1-RF_samps).*Bank_samps(:,tind);
    Conc_samps(:,t) = exp(-1/LT(end))*Conc_samps(:,t-1)+(1/ppt_to_tonnes(mol_ii))*emissions(:,t);
    tind = tind+1;
end
Conc_LB(:,scen_ii,mol_ii) = Conc(:,scen_ii,mol_ii)-prctile(Conc_samps,50)';

subplot(1,3,3)
scen_ii = 2;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.2,0.2,0.2]);
p1 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'Color',[0.2,0.2,0.2],'LineWidth',2);
hold on; 
scen_ii = 3;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.8500, 0.3250, 0.0980]);
hold on; 
p2 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
hold on; 
p3 = plot([y1:y3],Conc(:,1,mol_ii),'Color',[0, 0.4470, 0.7410],'LineWidth',2);
ylabel('CFC-113 [ppt]'); ylim([0,120]);

box on; 
legend([p1 p2 p3], 'Scen 1', 'Scen2','WMO 2018')
set(gca,'FontSize',12);

figure_width = 18; % in inches
figure_height = 6; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig5.pdf');
print(gcf, '-dpdf', str);


%% With inserts
FigHandle = figure(13); 
subplot(1,3,1)
mol_ii = 1;
scen_ii = 2;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.2,0.2,0.2]);
hold on; 
p1 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'Color',[0.2,0.2,0.2],'LineWidth',2)
scen_ii = 3;
hold on; 
p2 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2);
%p2 = boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0,0,0.7]);
hold on; 
scen_ii = 4;
p3 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2);
hold on; 
p4 = plot([y1:y3],Conc(:,1,mol_ii),'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
xlim([2010,2040]);
box on; set(gca,'FontSize',12);

subplot(1,3,2)
mol_ii = 2;
scen_ii = 2;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.2,0.2,0.2]);
hold on;
p1 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2,'Color',[0.2,0.2,0.2]);
hold on;
p2 = plot([y1:y3],Conc(:,1,mol_ii),'LineWidth',2,'Color',[0, 0.4470, 0.7410]);
xlim([2010,2040]);
box on; set(gca,'FontSize',12);

subplot(1,3,3)
mol_ii = 3;
scen_ii = 2;
boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.2,0.2,0.2]);
p1 = plot([y1:y3],Conc(:,scen_ii,mol_ii),'Color',[0.2,0.2,0.2],'LineWidth',2);
hold on; 
scen_ii = 3;
p2 = boundedline([y1:y3],Conc(:,scen_ii,mol_ii),[Conc_LB(:,scen_ii,mol_ii),Conc_UB(:,scen_ii,mol_ii)],'alpha','cmap',[0.8500, 0.3250, 0.0980]);
hold on; plot([y1:y3],Conc(:,scen_ii,mol_ii),'LineWidth',2,'Color',[0.8500, 0.3250, 0.0980]);
hold on; 
p3 = plot([y1:y3],Conc(:,1,mol_ii),'Color',[0, 0.4470, 0.7410],'LineWidth',2);
xlim([2010,2040]);
 box on; set(gca,'FontSize',12); 
 
 figure_width = 7; % in inches
figure_height = 1.8; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/Fig5_insert.pdf');
print(gcf, '-dpdf', str);
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Figures for Supplement %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
close all
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

%% Figure S1: Sparc Lifetimes 


load('SPARC_lifetimes.mat')
FigHandle = figure(1)
mol_names{1} = 'CFC-11';
mol_names{2} = 'CFC-12';
mol_names{3} = 'CFC-113';

Mod_names{1} =  'GSFC2D';   
Mod_names{2} =  'GEOSCCM';
Mod_names{3} =  'LMDZrepro';     
Mod_names{4} =  'SOCOL';
Mod_names{5} =  'ULAQ';
Mod_names{6} =  'UMUKCA';
Mod_names{7} =  'WACCM';

const_ii = [45, 100, 85];

year = 1960:2010;

for mol_ii = 1:3
    subplot(3,1,mol_ii)
    tmp = squeeze(LT_sparc(:,:,mol_ii));
    p2 = plot(year,tmp,'LineWidth',1.5);
    hold on; 
    p1 = plot(year,nanmean(tmp'),'k','LineWidth',3);
    title(mol_names{mol_ii}); ylabel('Lifetime [yrs]')
    p3 = plot(year,const_ii(mol_ii)*ones(size(year)), '--r','LineWidth',2);
    if mol_ii == 1
        legend([p1; p2],'MMM',Mod_names{:});
        ylim([40,120])
    end
    set(gca,'FontSize',11); 
    hold on; 

end
xlabel('Year')

figure_width = 7; % in inches
figure_height = 9; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig1.pdf');
print(gcf, '-dpdf', str);

%% Supplement Fig 2: Prior and posterior DE for MMM and constant LT

clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
load(CFC11_filename.meanLT);

var1 = DE_prior;
var3 = DE_prior(:,ResampleIndex);
FigHandle  = figure(2);
clear tmp
clear HBINS
for ii = 1:2
subplot(2,2,2*(ii-1)+1);
indx = randi(size(var1,2),100000,1);

tmp(1,:) = var1(20*ii+10,indx);
indx = randi(size(var3,2),100000,1);
tmp(2,:) = var3(20*ii+10,indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
HBINS(ii,:) = h.BinEdges;
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
str = num2str(y1+20*ii+10); title(str); xlabel('DE'); ylabel('Count');
legend('prior','posterior (SPARC MMM LT)'); set(gca, 'FontSize',12); 
end


load(CFC11_filename.LT45);

var1 = DE_prior;
var3 = DE_prior(:,ResampleIndex);
clear tmp
for ii = 1:2
subplot(2,2,2*ii);
indx = randi(size(var1,2),100000,1);
tmp(1,:) = var1(20*ii+10,indx);
indx = randi(size(var3,2),100000,1);
tmp(2,:) = var3(20*ii+10,indx);
histogram(tmp(1,:)',HBINS(ii,:),'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',HBINS(ii,:),'facecolor',[1 0 0],'facealpha',.5);
str = num2str(y1+20*ii+10); title(str); xlabel('DE'); ylabel('Count');
legend('prior','posterior (LT = 45 yrs)'); set(gca, 'FontSize',12); 
end


figure_width = 18; % in inches
figure_height = 10; % in inches
screen_ppi = 72;

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig2.pdf');
print(gcf, '-dpdf', str);
%% Supplement Fig 3: Prior and posterior RF for MMM and constant LT
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
load(CFC11_filename.meanLT);

var1 = RF_prior;
var3 = RF_prior(:,ResampleIndex);
FigHandle  = figure(3);

for ii = 1:2
subplot(2,2,2*(ii-1)+1);
indx = randi(size(var1,2),100000,1);

tmp(1,:) = var1(20*ii+10,indx);
indx = randi(size(var3,2),100000,1);
tmp(2,:) = var3(20*ii+10,indx);
h = histogram(tmp(1,:)',40,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
HBINS(ii,:) = h.BinEdges;
str = num2str(y1+20*ii+10); title(str); xlabel('RF'); ylabel('Count');
legend('prior','posterior (SPARC MMM LT)'); set(gca, 'FontSize',12); 
end

load(CFC11_filename.LT45);

var1 = RF_prior;
var3 = RF_prior(:,ResampleIndex);

for ii = 1:2
subplot(2,2,2*ii);
indx = randi(size(var1,2),100000,1);
tmp(1,:) = var1(20*ii+10,indx);
indx = randi(size(var3,2),100000,1);
tmp(2,:) = var3(20*ii+10,indx);
histogram(tmp(1,:)',HBINS(ii,:),'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',HBINS(ii,:),'facecolor',[1 0 0],'facealpha',.5);
str = num2str(y1+20*ii+10); title(str); xlabel('RF'); ylabel('Count');
legend('prior','LT = 45 yrs'); set(gca, 'FontSize',12); 
end


figure_width = 18; % in inches
figure_height = 10; % in inches
screen_ppi = 72;

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);


str = strcat(FolderName,'/SI_Fig3.pdf');
print(gcf, '-dpdf', str);

%% Supplement Fig 4,5,6, and 7: Joint priors and posteriors between DE,RF, prod and Bank
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
load(CFC11_filename.meanLT);
var1 = DE_prior(:,ResampleIndex);
var2 = RF_prior(:,ResampleIndex);
var3 = Prod_prior(ResampleIndex,:)';
var4 = 0.001*Bank(ResampleIndex,:)';
var5 = 0.001*Emiss(ResampleIndex,:)';
var1prior = DE_prior;
var2prior = RF_prior;
var3prior = Prod_prior';
var4prior = 0.001*Bank';
var5prior = 0.001*Emiss';




FigHandle  = figure(4);
clrmp = cbrewer('div','RdGy',80);
clrmp = clrmp(1:30,:);
clrmp = flipud(clrmp);
ystr = 'Bank [Gg]';
for ii = 1:3
    subplot(2,3,ii);
    xstr = 'DE';
    h = histogram2(var1prior(10*ii+20,:)',var4prior(10*ii+20,:)',[80,80],'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint Bank and DE Prior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12)
    
    subplot(2,3,ii+3);
    xstr = 'DE';
    histogram2(var1(10*ii+20,:)',var4(10*ii+20,:)',h.XBinEdges,h.YBinEdges,'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint Bank and Posterior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12)

end

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72;

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig4.pdf');
print(gcf, '-dpdf', str);



FigHandle  = figure(5);
ystr = 'Bank [Gg]';
for ii = 1:3
    subplot(2,3,ii);
    xstr = 'RF';
    h = histogram2(var2prior(10*ii+20,:)',var4prior(10*ii+20,:)',[80,80],'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint Bank and RF Prior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12); 
    %hold on; plot(h.XBinEdges, repmat(0.001*ObsDerivedEmiss(10*ii+20),1,81)./h.XBinEdges,'k')
    
    subplot(2,3,ii+3);
    xstr = 'RF';
    histogram2(var2(10*ii+20,:)',var4(10*ii+20,:)',h.XBinEdges,h.YBinEdges,'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint Bank and RF Posterior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12); 
    %hold on; plot(h.XBinEdges, repmat(0.001*ObsDerivedEmiss(10*ii+20),1,81)./h.XBinEdges,'k','LineWidth',2)
end

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig5.pdf');
print(gcf, '-dpdf', str);


FigHandle  = figure(6);
ystr = 'Bank [Gg]';
for ii = 1:3
    subplot(2,3,ii);
    xstr = 'Production [Gg]';
    h = histogram2(0.001*var3prior(10*ii+20,:)',var4prior(10*ii+20,:)',[80,80],'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint Bank and Prod Prior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12)
    
    subplot(2,3,ii+3);
    xstr = 'Production [Gg]';
    histogram2(0.001*var3(10*ii+20,:)',var4(10*ii+20,:)',h.XBinEdges,h.YBinEdges,'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint Bank and Prod Posterior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12)

end

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig6.pdf');
print(gcf, '-dpdf', str);


FigHandle  = figure(7);
ystr = 'DE';
for ii = 1:3
    subplot(2,3,ii);
    xstr = 'RF';
    h = histogram2(var2prior(10*ii+20,:)',var1prior(10*ii+20,:)',[80,80],'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint RF and DE Prior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12)
    
    subplot(2,3,ii+3);
    histogram2(var2(10*ii+20,:)',var1(10*ii+20,:)',h.XBinEdges,h.YBinEdges,'DisplayStyle','tile','Normalization','probability','LineStyle','none'); 
    xlabel(xstr); ylabel(ystr); shading flat; grid off; 
    str = num2str(y1+10*ii+20); str = strcat('Joint RF and DE Posterior, ',str); title(str); colorbar;
    colormap(clrmp); caxis([0,0.015]);
    xlim([h.XBinEdges(1)-h.BinWidth(1),h.XBinEdges(end)+h.BinWidth(1)]);
    ylim([h.YBinEdges(1)-h.BinWidth(2),h.YBinEdges(end)+h.BinWidth(2)]);
    set(gca, 'FontSize',12)

end

figure_width = 18; % in inches
figure_height = 9; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig7.pdf');
print(gcf, '-dpdf', str);

%% Supplement Fig 8: Prior and posterior Production for MMM and constant LT

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
load(CFC11_filename.meanLT);

var1 = Prod_prior';
var3 = Prod_prior(ResampleIndex,:)';
FigHandle  = figure(8);
clear tmp
clear HBINS
HBINS(1,:) = 340000:10000:650000;
HBINS(2,:) = 21300:275:30000;
for ii = 1:2
subplot(2,2,2*(ii-1)+1);
indx = randi(size(var1,2),100000,1);

tmp(1,:) = var1(20*ii+10,indx);
indx = randi(size(var3,2),100000,1);
tmp(2,:) = var3(20*ii+10,indx);
histogram(tmp(1,:)',HBINS(ii,:),'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',HBINS(ii,:),'facecolor',[1 0 0],'facealpha',.5);
str = num2str(y1+20*ii+10); title(str); xlabel('Production'); ylabel('Count');
legend('prior','posterior (SPARC MMM LT)'); set(gca, 'FontSize',12); 
end

load(CFC11_filename.LT45);
clear tmp
var1 = Prod_prior';
var3 = Prod_prior(ResampleIndex,:)';

for ii = 1:2
subplot(2,2,2*ii);
indx = randi(size(var1,2),100000,1);
tmp(1,:) = var1(20*ii+10,indx);
indx = randi(size(var3,2),100000,1);
tmp(2,:) = var3(20*ii+10,indx);
histogram(tmp(1,:)',HBINS(ii,:),'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',HBINS(ii,:),'facecolor',[1 0 0],'facealpha',.5);
str = num2str(y1+20*ii+10); title(str); xlabel('Production'); ylabel('Count');
legend('prior','posterior (fixed LT)'); set(gca, 'FontSize',12); 
end


figure_width = 18; % in inches
figure_height = 10; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);


str = strcat(FolderName,'/SI_Fig8.pdf');
print(gcf, '-dpdf', str);
%% SI Fig 9: Prior and posterior uncertainties for each molecule
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
FigHandle = figure(9)
subplot(2,3,1)
load(CFC11_filename.meanLT);
var1 = emiss_sigma_prior';
var3 = emiss_sigma_prior(ResampleIndex,:)';
clear tmp
indx = randi(size(var1,2),10000,1);
tmp(1,:) = 0.001*var1(10,indx);
indx = randi(size(var3,2),10000,1);
tmp(2,:) = 0.001*var3(10,indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-11'); xlabel('\sigma \times UB (emissions [Gg])'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(2,3,4)
tmp(1,:) = 1000*(1/mean(ObsDerivedEmiss))*tmp(1,:);
tmp(2,:) = 1000*(1/mean(ObsDerivedEmiss))*tmp(2,:);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-11'); xlabel('relative uncertainty'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(2,3,2)
load(CFC12_filename.meanLT);
var1 = emiss_sigma_prior';
var3 = emiss_sigma_prior(ResampleIndex,:)';
clear tmp
indx = randi(size(var1,2),10000,1);
tmp(1,:) = 0.001*var1(10,indx);
indx = randi(size(var3,2),10000,1);
tmp(2,:) = 0.001*var3(10,indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-12'); xlabel('\sigma \times UB  (emissions [Gg])'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(2,3,5)
tmp(1,:) = 1000*(1/mean(ObsDerivedEmiss))*tmp(1,:);
tmp(2,:) = 1000*(1/mean(ObsDerivedEmiss))*tmp(2,:);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-12'); xlabel('relative uncertainty'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(2,3,3)
load(CFC113_filename.meanLT);
var1 = emiss_sigma_prior';
var3 = emiss_sigma_prior(ResampleIndex,:)';
clear tmp
indx = randi(size(var1,2),10000,1);
tmp(1,:) = 0.001*var1(10,indx);
indx = randi(size(var3,2),10000,1);
tmp(2,:) = 0.001*var3(10,indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-113'); xlabel('\sigma \times UB (emissions [Gg])'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(2,3,6)
tmp(1,:) = 1000*(1/mean(ObsDerivedEmiss))*tmp(1,:);
tmp(2,:) = 1000*(1/mean(ObsDerivedEmiss))*tmp(2,:);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-12'); xlabel('relative uncertainty'); ylabel('Count');
set(gca, 'FontSize',12)

figure_width = 22; % in inches
figure_height = 12; % in inches
screen_ppi = 72;

screen_figure_width = round(figure_width*screen_ppi); % in pixels

screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig9.pdf');
print(gcf, '-dpdf', str);


%% Supplement Fig 10: Residuals
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
FigHandle = figure(10)
load(CFC11_filename.meanLT);

yrindx = [1980,1990,2000,2010];
for ii = 1:4
    subplot(2,4,ii); 
    ytmp1 = find(years == yrindx(ii));
    hist(Emiss(ResampleIndex,ytmp1)-ObsDerivedEmiss(ytmp1),30); 
    title(num2str(yrindx(ii))); xlabel('residuals [tonnes CFC-11]');
    set(gca, 'FontSize',12)
    subplot(2,4,ii+4);    
    qqplot(Emiss(ResampleIndex,ytmp1)-ObsDerivedEmiss(ytmp1)); 
    title(num2str(yrindx(ii))); box on; set(gca, 'FontSize',12)
end

figure_width = 22; % in inches
figure_height = 6; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig10.pdf');
print(gcf, '-dpdf', str);

%% Updated Number Fig S11 Emissions posterior with unexpected production scenario

clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;


SCEN_FILE{1} = CFC11_filename.fugitive;
Mole_Name{1} = 'CFC-11';
Input_FileName{1} = 'CFC11';
LT_const(1) = 62.9;
ppt_to_tonnes(1) =22602.38457; % Conversion of ppt to tonnes of CFC-11

FigHandle = figure(11);
scen_ii = 1
load(SCEN_FILE{scen_ii})
MED = prctile(Emiss(ResampleIndex,:),50);
UB = (prctile(Emiss(ResampleIndex,:),97.5)-MED);
LB = (MED-prctile(Emiss(ResampleIndex,:),2.5));
p1 = boundedline([y1:y2],0.001*MED',[0.001*LB',0.001*UB'],'alpha','cmap',[0.3,0.3,0.3]);
UBmax = (prctile(Emiss(ResampleIndex,:),99.5));
LBmin = (prctile(Emiss(ResampleIndex,:),0.5));
hold on; plot([y1:y2],0.001*UBmax,'--k');
hold on; plot([y1:y2],0.001*LBmin,'--k');
hold on; 
p2 = plot([y1:y2],0.001*ObsDerivedEmiss,'r','LineWidth',2);
str = strcat(Mole_Name{scen_ii}, ' Emissions');
title(str);
str = strcat(Mole_Name{scen_ii}, ' emissions [Gg yr^{-1}]');
box on; ylabel(str); xlabel('year'); ylim([0,0.001*1.7*max(ObsDerivedEmiss)])
set(gca,'FontSize',12)
xlim([1955,2016])

str = strcat(Input_FileName{scen_ii},'/Input/wmo2018.mat');
load(str);
ytmp1 = find(wmo_yr == y1);
ytmp2 = find(wmo_yr == y2);
conc_obs = wmo_conc(ytmp1:ytmp2+1);
Emiss_constantLT = ppt_to_tonnes(scen_ii)*(conc_obs(2:end)'-conc_obs(1:end-1)'.*exp(-1./LT_const(scen_ii)'));
hold on; 
p3 = plot([y1:y2],0.001*Emiss_constantLT,'g','LineWidth',2);
str = strcat('LT =',' ',num2str(LT_const(scen_ii)),' yrs');
legend([p1 p2 p3],{'bayesian posterior','LT = MMM', str})


figure_width = 8; % in inches
figure_height = 5; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig11.pdf');
print(gcf, '-dpdf', str);

%% Sensitivity of Banks to RF, DE and Prod
% To create this figure, first rerun MainScript_multi with the following
% specifications: 
% CFC-11, Make Prod_mu 10% higher (name run CFC11_prod10pct_meanSparcLTrun.mat
% CFC-11, First, rerun the joint RF and DE priors, with RF values 10%
% higher for each equipment type.  then run MainScript_multi.m and save
% output as CFC11_RF10pct_meanSparcLTrun.mat'
% CFC-11, First, rerun the joint RF and DE priors, with DE values 10%
% higher for each equipment type.  Then run MainScript_multi.m and save
% output as CFC11_RF10pct_meanSparcLTrun.mat'
% CFC-11, First, rerun the joint RF and DE priors, with RF standard deviation of in
% closed cell foams being twice as large. THen run MainScript_multi.m and
% save output as 'CFC11/Output/CFC11_RFhyper_meanSparcLTrun.mat'. 


close all
clear all
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

CFC11_filename.meanLT = strcat(HomeDir,'CFC11/Output/CFC11_meanSparcLT.mat');
CFC11_filename.Prod10 = strcat(HomeDir,'CFC11/Output/CFC11_prod10pct_meanSparcLTrun.mat';
CFC11_filename.RF10 = strcat(HomeDir,'CFC11/Output/CFC11_RF10pct_meanSparcLTrun.mat';
CFC11_filename.DE10 = strcat(HomeDir,'CFC11/Output/CFC11_DE10pct_meanSparcLTrun.mat';
CFC11_filename.RFunc = strcat(HomeDir,'CFC11/Output/CFC11_RFhyper_meanSparcLTrun.mat';


load(CFC11_filename.meanLT,'Bank','Banks_TD','ResampleIndex','Emiss','ObsDerivedEmiss')
Bank_BPE = Bank(ResampleIndex,:);
MED_bc = prctile(Bank_BPE,50);
UB_bc = (prctile(Bank_BPE,97.5)-MED_bc)';
LB_bc = (MED_bc-prctile(Bank_BPE,2.5))';
Bank_meanLT = Banks_TD;

Emiss_BPE = Emiss(ResampleIndex,:);
Emiss_MED_bc = prctile(Emiss_BPE,50);
Emiss_UB_bc = (prctile(Emiss_BPE,97.5)-Emiss_MED_bc)';
Emiss_LB_bc = (Emiss_MED_bc-prctile(Emiss_BPE,2.5))';


FigHandle = figure; 

subplot(2,2,1)
p1 = boundedline([y1:y2],0.001*MED_bc',0.001*[LB_bc,UB_bc],'alpha','cmap',[0.7,0,0]); hold on; 
p2 = plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2);

load(CFC11_filename.Prod10,'Bank','Banks_TD','ResampleIndex','Emiss','ObsDerivedEmiss')
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50); UB = (prctile(Bank_BPE,97.5)-MED)'; LB = (MED-prctile(Bank_BPE,2.5))';
hold on; 
p3 = plot([y1:y2],0.001*Banks_TD,'--b','LineWidth',2); hold on; 
p4 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]);
xlim([y1,y2]); xlabel('year'); ylabel('CFC-11 Bank [Gg]');
lgd = legend([p1 p4 p2 p3],'base case','prod+10%', 'base case TD','prod+10% TD'); set(gca,'FontSize',13); 
lgd.Location = 'northwest'; box on; title('Prod + 10%')

subplot(2,2,2)

p1 = boundedline([y1:y2],0.001*Emiss_MED_bc',0.001*[Emiss_LB_bc,Emiss_UB_bc],'alpha','cmap',[0.7,0,0]); hold on;

Emiss_BPE = Emiss(ResampleIndex,:);
MED = prctile(Emiss_BPE,50);
UB = (prctile(Emiss_BPE,97.5)-MED)';
LB = (MED-prctile(Emiss_BPE,2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]); hold on;
p3 = plot([y1:y2],0.001*ObsDerivedEmiss,'--k','LineWidth',2);
title('posterior emissions space'); box on; ylabel('Emissions [Gg/yr]');
set(gca,'FontSize',13); xlim([1955,2016]); xlabel('year');
lgd = legend([p1 p2 p3],'base case','prod+10%','obs-derived'); set(gca,'FontSize',13); 
lgd.Location = 'northwest'; box on; title('Prod + 10%')


subplot(2,2,3)
p1 = boundedline([y1:y2],0.001*MED_bc',0.001*[LB_bc,UB_bc],'alpha','cmap',[0.7,0,0]); hold on;

load(CFC11_filename.DE10,'Bank','Banks_TD','ResampleIndex','Emiss','ObsDerivedEmiss')
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50); UB = (prctile(Bank_BPE,97.5)-MED)'; LB = (MED-prctile(Bank_BPE,2.5))';
p3 = plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2); hold on; 
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]);
xlim([y1,y2]); xlabel('year'); ylabel('CFC-11 Bank [Gg]'); title('DE + 10%')
lgd = legend([p1 p2 p3],'base case','DE+10%','top-down'); box on;
lgd.Location = 'northwest'; set(gca,'FontSize',13); 

subplot(2,2,4)
p1 = boundedline([y1:y2],0.001*Emiss_MED_bc',0.001*[Emiss_LB_bc,Emiss_UB_bc],'alpha','cmap',[0.7,0,0]); hold on; 
Emiss_BPE = Emiss(ResampleIndex,:);
MED = prctile(Emiss_BPE,50);
UB = (prctile(Emiss_BPE,97.5)-MED)';
LB = (MED-prctile(Emiss_BPE,2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]); hold on;
p3 = plot([y1:y2],0.001*ObsDerivedEmiss,'--k','LineWidth',2);
title('posterior emissions space'); box on; ylabel('Emissions [Gg/yr]');
set(gca,'FontSize',13); xlim([1955,2016]); xlabel('year'); title('DE + 10%');
lgd = legend([p1 p2 p3],'base case','DE+10%','obs-derived'); set(gca,'FontSize',13); 
lgd.Location = 'northwest'; box on; 


figure_width = 16; % in inches
figure_height = 8; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig12TOP.pdf');
print(gcf, '-dpdf', str);

FigHandle = figure; 
subplot(2,2,1)
p1 = boundedline([y1:y2],0.001*MED_bc',0.001*[LB_bc,UB_bc],'alpha','cmap',[0.7,0,0]);
hold on; plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2)

load(CFC11_filename.RF10,'Bank','Banks_TD','ResampleIndex','Emiss','ObsDerivedEmiss')
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50); UB = (prctile(Bank_BPE,97.5)-MED)'; LB = (MED-prctile(Bank_BPE,2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]);
xlim([y1,y2]); xlabel('year'); ylabel('CFC-11 Bank [Gg]'); title('RF + 10%')
lgd = legend([p1 p2],'base case','RF+10%');
lgd.Location = 'northwest'; box on; set(gca,'FontSize',13); 


subplot(2,2,2)
p1 = boundedline([y1:y2],0.001*Emiss_MED_bc',0.001*[Emiss_LB_bc,Emiss_UB_bc],'alpha','cmap',[0.7,0,0]);
Emiss_BPE = Emiss(ResampleIndex,:);
MED = prctile(Emiss_BPE,50);
UB = (prctile(Emiss_BPE,97.5)-MED)';
LB = (MED-prctile(Emiss_BPE,2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]); hold on;
hold on; 
p3 = plot([y1:y2],0.001*ObsDerivedEmiss,'--k','LineWidth',2);
title('posterior emissions space'); box on; ylabel('Emissions [Gg/yr]')
set(gca,'FontSize',13); xlim([1955,2016]);xlabel('year');
lgd = legend([p1 p2 p3],'base case','RF+10%','obs-derived'); set(gca,'FontSize',13); 
lgd.Location = 'northwest'; box on; title('RF + 10%')


subplot(2,2,3)
p1 = boundedline([y1:y2],0.001*MED_bc',0.001*[LB_bc,UB_bc],'alpha','cmap',[0.7,0,0]);
hold on; plot([y1:y2],0.001*Banks_TD,'--r','LineWidth',2)

load(CFC11_filename.RFunc,'Bank','Banks_TD','ResampleIndex','Emiss','ObsDerivedEmiss')
Bank_BPE = Bank(ResampleIndex,:);
MED = prctile(Bank_BPE,50); UB = (prctile(Bank_BPE,97.5)-MED)'; LB = (MED-prctile(Bank_BPE,2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]);
xlim([y1,y2]); xlabel('year'); ylabel('CFC-11 Bank [Gg]'); title('sensitivity to RF uncertainty')
lgd = legend([p1 p2],'base case','RFunc x2');
lgd.Location = 'northwest'; box on; set(gca,'FontSize',13); 

subplot(2,2,4)
p1 = boundedline([y1:y2],0.001*Emiss_MED_bc',0.001*[Emiss_LB_bc,Emiss_UB_bc],'alpha','cmap',[0.7,0,0]);
Emiss_BPE = Emiss(ResampleIndex,:);
MED = prctile(Emiss_BPE,50);
UB = (prctile(Emiss_BPE,97.5)-MED)';
LB = (MED-prctile(Emiss_BPE,2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0,0,0.7]); hold on;
hold on; 
p3 = plot([y1:y2],0.001*ObsDerivedEmiss,'--k','LineWidth',2);
title('RF unc x2'); box on; ylabel('Emissions [Gg/yr]');
set(gca,'FontSize',13); xlim([1955,2016]);xlabel('year');
lgd = legend([p1 p2 p3],'base case','RF unc x 2%','obs-derived'); set(gca,'FontSize',13); 
lgd.Location = 'northwest'; box on; title('RF unc x 2')


figure_width = 16; % in inches
figure_height = 8; % in inches
screen_ppi = 72;

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

FolderName = 'April3';
str = strcat(FolderName,'/SI_Fig12BOTTOM.pdf');
print(gcf, '-dpdf', str);

%% Supplement Fig 14: Timeseries of priors and posteriors for CFC-11
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

load(CFC11_filename.meanLT)

FigHandle = figure(14);
subplot(3,1,1)
MED = prctile(DE_prior',50);
UB = (prctile(DE_prior',97.5)-MED)';
LB = (MED-prctile(DE_prior',2.5))';
p1 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
hold on; 
MED = prctile(DE_prior(:,ResampleIndex)',50);
UB = (prctile(DE_prior(:,ResampleIndex)',97.5)-MED)';
LB = (MED-prctile(DE_prior(:,ResampleIndex)',2.5))';
p2 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('DE'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-11 Direct Emissions')
legend([p1 p2],'Prior DE','Posterior DE');

subplot(3,1,2)
MED = prctile(RF_prior',50);
UB = (prctile(RF_prior',97.5)-MED)';
LB = (MED-prctile(RF_prior',2.5))';
p1 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
hold on; 
MED = prctile(RF_prior(:,ResampleIndex)',50);
UB = (prctile(RF_prior(:,ResampleIndex)',97.5)-MED)';
LB = (MED-prctile(RF_prior(:,ResampleIndex)',2.5))';
p2 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('RF'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-11 Release Fraction')
legend([p1 p2],'Prior RF','Posterior RF');


subplot(3,1,3)
MED = prctile(Prod_prior,50);
UB = (prctile(Prod_prior,97.5)-MED)';
LB = (MED-prctile(Prod_prior,2.5))';
p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
hold on; 
MED = prctile(Prod_prior(ResampleIndex,:),50);
UB = (prctile(Prod_prior(ResampleIndex,:),97.5)-MED)';
LB = (MED-prctile(Prod_prior(ResampleIndex,:),2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('Production [Gg]'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-11 Production')
legend([p1 p2],'Prior Production','Posterior Production');

figure_width = 8; % in inches
figure_height = 12; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig14.pdf');
print(gcf, '-dpdf', str);

%% Supplement Fig 15: Timeseries of priors and posteriors for CFC-12
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

load(CFC12_filename.meanLT)

FigHandle = figure(15);

subplot(3,1,1)
MED = prctile(DE_prior',50);
UB = (prctile(DE_prior',97.5)-MED)';
LB = (MED-prctile(DE_prior',2.5))';
p1 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
hold on; 
MED = prctile(DE_prior(:,ResampleIndex)',50);
UB = (prctile(DE_prior(:,ResampleIndex)',97.5)-MED)';
LB = (MED-prctile(DE_prior(:,ResampleIndex)',2.5))';
p2 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('DE'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-12 Direct Emissions')
legend([p1 p2],'Prior DE','Posterior DE');

subplot(3,1,2)
MED = prctile(RF_prior',50);
UB = (prctile(RF_prior',97.5)-MED)';
LB = (MED-prctile(RF_prior',2.5))';
p1 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
hold on; 
MED = prctile(RF_prior(:,ResampleIndex)',50);
UB = (prctile(RF_prior(:,ResampleIndex)',97.5)-MED)';
LB = (MED-prctile(RF_prior(:,ResampleIndex)',2.5))';
p2 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('RF'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-12 Release Fraction')
legend([p1 p2],'Prior RF','Posterior RF');

subplot(3,1,3)
MED = prctile(Prod_prior,50);
UB = (prctile(Prod_prior,97.5)-MED)';
LB = (MED-prctile(Prod_prior,2.5))';
p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
hold on; 
MED = prctile(Prod_prior(ResampleIndex,:),50);
UB = (prctile(Prod_prior(ResampleIndex,:),97.5)-MED)';
LB = (MED-prctile(Prod_prior(ResampleIndex,:),2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('Production [Gg]'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-12 Production')
legend([p1 p2],'Prior Production','Posterior Production');

figure_width = 8; % in inches
figure_height = 12; % in inches
screen_ppi = 72;

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig15.pdf');
print(gcf, '-dpdf', str);

%% Supplement Fig S16: Timeseries of priors and posteriors for CFC-113
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;

load(CFC113_filename.meanLT)

FigHandle = figure(16);
subplot(3,1,1)
MED = prctile(DE_prior',50);
UB = (prctile(DE_prior',97.5)-MED)';
LB = (MED-prctile(DE_prior',2.5))';
p1 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
%boundedline([y1:y2],0.001*mean(Posterior_emiss_samps,2),0.001*1.96*[std(Posterior_emiss_samps')',std(Posterior_emiss_samps')'],'alpha','cmap',[0,0,1]);
hold on; 
MED = prctile(DE_prior(:,ResampleIndex)',50);
UB = (prctile(DE_prior(:,ResampleIndex)',97.5)-MED)';
LB = (MED-prctile(DE_prior(:,ResampleIndex)',2.5))';
p2 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('DE'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-113 Direct Emissions')
legend([p1 p2],'Prior DE','Posterior DE');

subplot(3,1,2)
MED = prctile(RF_prior',50);
UB = (prctile(RF_prior',97.5)-MED)';
LB = (MED-prctile(RF_prior',2.5))';
p1 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
%boundedline([y1:y2],0.001*mean(Posterior_emiss_samps,2),0.001*1.96*[std(Posterior_emiss_samps')',std(Posterior_emiss_samps')'],'alpha','cmap',[0,0,1]);
hold on; 
MED = prctile(RF_prior(:,ResampleIndex)',50);
UB = (prctile(RF_prior(:,ResampleIndex)',97.5)-MED)';
LB = (MED-prctile(RF_prior(:,ResampleIndex)',2.5))';
p2 = boundedline([y1:y2],MED',[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('RF'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-113 Release Fraction')
legend([p1 p2],'Prior RF','Posterior RF');

subplot(3,1,3)
MED = prctile(Prod_prior,50);
UB = (prctile(Prod_prior,97.5)-MED)';
LB = (MED-prctile(Prod_prior,2.5))';
p1 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[0.3,0.3,0.3]);
%boundedline([y1:y2],0.001*mean(Posterior_emiss_samps,2),0.001*1.96*[std(Posterior_emiss_samps')',std(Posterior_emiss_samps')'],'alpha','cmap',[0,0,1]);
hold on; 
MED = prctile(Prod_prior(ResampleIndex,:),50);
UB = (prctile(Prod_prior(ResampleIndex,:),97.5)-MED)';
LB = (MED-prctile(Prod_prior(ResampleIndex,:),2.5))';
p2 = boundedline([y1:y2],0.001*MED',0.001*[LB,UB],'alpha','cmap',[1,0,0]);
ylabel('Production [Gg]'); set(gca, 'FontSize',12); xlabel('Year');
xlim([y1,y2]); box on; title('CFC-113 Production')
legend([p1 p2],'Prior Production','Posterior Production');

figure_width = 8; % in inches
figure_height = 12; % in inches
screen_ppi = 72; 

screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig16.pdf');
print(gcf, '-dpdf', str);

%% SI Figure 18: Production Correlation
clearvars -except CFC11_filename CFC12_filename CFC113_filename FolderName

y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
FigHandle = figure(18)
subplot(1,3,1)
load(CFC11_filename.meanLT);
var1 = squeeze(rho1);
var3 = squeeze(rho1(1,1,ResampleIndex));
clear tmp
indx = randi(size(var1,1),10000,1);
tmp(1,:) = var1(indx);
indx = randi(size(var3,1),10000,1);
tmp(2,:) = var3(indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-11'); xlabel('\rho'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(1,3,2)
load(CFC12_filename.meanLT);
var1 = squeeze(rho1);
var3 = squeeze(rho1(1,1,ResampleIndex));
clear tmp
indx = randi(size(var1,1),10000,1);
tmp(1,:) = var1(indx);
indx = randi(size(var3,1),10000,1);
tmp(2,:) = var3(indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-12'); xlabel('\rho'); ylabel('Count');
set(gca, 'FontSize',12)

subplot(1,3,3)
load(CFC113_filename.meanLT);
var1 = squeeze(rho1);
var3 = squeeze(rho1(1,1,ResampleIndex));
clear tmp
indx = randi(size(var1,1),10000,1);
tmp(1,:) = var1(indx);
indx = randi(size(var3,1),10000,1);
tmp(2,:) = var3(indx);
h = histogram(tmp(1,:)',25,'facecolor',[0 0 1]);
hold on; 
histogram(tmp(2,:)',h.BinEdges,'facecolor',[1 0 0],'facealpha',.5);
title('CFC-113'); xlabel('\rho'); ylabel('Count');
set(gca, 'FontSize',12)


figure_width = 22; % in inches
figure_height = 6; % in inches
%font_size = 6; % in points (relative to figure size in inches)
export_ppi = 300; % % number of pixels per inch *in exported file*; only relevant if exporting as pixel-based figure (such as png, tiff)

% WE ALSO NEED
screen_ppi = 72; % inherent property of Matlab. CHANGING THIS WILL MAKE FIGURE LOOK INCONSISTENT BETWEEN SCREEN AND EXPORT!

% DERIVED PROPERTIES (cannot be changed; for info only)
screen_figure_width = round(figure_width*screen_ppi); % in pixels
screen_figure_height = round(figure_height*screen_ppi); % in pixels
set(FigHandle, 'Position', [100, 100, round(figure_width*screen_ppi), round(figure_height*screen_ppi)]);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperSize', [figure_width figure_height]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 figure_width figure_height]);

str = strcat(FolderName,'/SI_Fig18.pdf');
print(gcf, '-dpdf', str);

