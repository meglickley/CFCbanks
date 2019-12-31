% BPE Banks model
% Last udpated Dec 1, 2019 by Megan Lickley.
% This script runs the BPE model 
% First run the estimating_RFandDE.m file. 
% Edit the HomeDir path
% Specify which lifetime scenario and which gas scenario to run
% If replicating all figures, then this script must be run for the
% following configurations:
        %  CFC-11 with mean SPARC lifetime scenario
        %  CFC-11 with constant lifetime scenario
        %  CFC-11 with constant lifetime, 2nd scenario (Constant_LT_scen2)
        %  CFC-11 with mean SPARC lifetime and Fugitive scenario
        %  CFC-12 with mean SPARC lifetime scenario
        %  CFC-12 with constant lifetime scenario
        %  CFC-113 with mean SPARC lifetime scenario
        %  CFC-113 with constant lifetime scenario
% Run this script for all scenarios then run Make_Figures.m 

close all
clear all
y1 = 1955; % Year Start Date
y2 = 2016; % Year End Date
years = y1:y2;
Nyears = length(years);
N = 5*10^5;
Nresamps = 10^5;

%% Change HomeDirectory name
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Code';

%% Specify Scenarios - change this depending on scenario and gas

% Specify which Life Time Scenario
meanSPARC = 1; % Runs the model with mean SPARC life time
Constant_LT = 0; % Assumes a constant life time throughout the simulation period
Constant_LT_scen2 = 0;

% Specify which gas 
CFC11 = 1;% Specify molecule
CFC12 = 0;
CFC113 = 0;

% Turn fugitive emissions on only for CFC11 for one case.  
fugitive_emissions = 1; 
%% Defining input parameters

if CFC11
    ppt_to_tonnes = 22602.38457; % Conversion of ppt to tonnes of CFC-11
    InputFileName = 'CFC11';
    Sparc_ID = 1;
    BankSize1 = 5893.9; % From WMO 2002
    const_LT_val = 45;%62.88; % Life time value - can change this to 45
    if Constant_LT_scen2
        const_LT_val = 62.88;
    end
elseif CFC12
    ppt_to_tonnes =19895.36; % Conversion of ppt to tonnes of CFC-12
    InputFileName = 'CFC12';
    Sparc_ID = 2;
    BankSize1 = 65198; % From WMO 2002
    const_LT_val = 100; % Constant lifetime value
elseif CFC113
    ppt_to_tonnes =30834.68811; % Conversion of ppt to tonnes of CFC-12
    InputFileName = 'CFC113';
    Sparc_ID = 3;
    BankSize1 = 1872; % From WMO 2002
    const_LT_val = 85; % Constant lifetime value
end

OutputFolderName = strcat(InputFileName,'/Output/');
%% Defining atmospherice life times
year = 1950:2099;
load('SPARC_lifetimes.mat')
tmp = squeeze(LT_sparc(:,:,Sparc_ID)); % Looking just at CFC-12
SLT = nan(150,size(tmp,2));
SLT(1:10,:) = repmat(tmp(1,:),10,1);
SLT(11:61,:) = tmp;
indx = ~isnan(tmp); % If data ends early, fill in the remainder with the last time period's
% entry
for modi = 1:size(SLT,2)
    SLT(sum(indx(:,modi))+1:end,modi) = SLT(sum(indx(:,modi)),modi);
end

if meanSPARC  
    LT = mean(SLT')';
    OutputFileName = strcat(InputFileName,'_meanSparcLT.mat');
else
    LT = const_LT_val*ones(150,1); 
    OutputFileName = strcat(InputFileName,'_fixedLT.mat');
    if Constant_LT_scen2
       OutputFileName = strcat(InputFileName,'_fixedLT62p9.mat');
    end
end

LT_min = min(SLT')';
LT_max = max(SLT')';
ytmp1 = find(year == y1); % Only keep years in desired range
ytmp2 = find(year == y2); 

LT = LT(ytmp1-1:ytmp2); 
LT_min = LT_min(ytmp1-1:ytmp2); 
LT_max = LT_max(ytmp1-1:ytmp2); 
%% Defining Priors for RF and DE - these come from the bottom-up accounting (see estimating_RFandDE.m)
tmp = strcat(InputFileName,'/Input');
cd(tmp)

if fugitive_emissions
    load('Fugitive_InferredRFandDE.mat'); 
else
    load('InferredRFandDE.mat');
end
ytmp1 = find(yr==y1);
ytmp2 = find(yr==y2);

RF_sim = RF(ytmp1:ytmp2,:);
DE_sim = DE(ytmp1:ytmp2,:);
RF_sim(RF_sim<0) = 0; % RF must be greater than zero. 

%% Define Priors for Production
Prod_mu = NaN(size(years)); % Production is in Tonnes
load('WMO2002.mat') % Use until 1988 - WMO data is the adjusted AFEAS data
ytmp1 = find(year==y1); 
ytmp2 = find(year==1988);
Prod_mu(1:1989-y1) = Production(ytmp1:ytmp2);

load('a5_na5.mat') % Use Article 5 and non A5 prod total from 1989 onwards
ytmp1 = find(a5_na5(:,1) == 1989);
ytmp2 = find(years == a5_na5(end,1));
Prod_mu(1989-y1+1:ytmp2) = sum(a5_na5(ytmp1:end,[2:3]),2);
Prod_mu(ytmp2:end) = sum(a5_na5(end,[2:3]),2);

% Parameters used to define the production prior for the simulation model
Prod_LL = Prod_mu;
Prod_UL = 1.1*Prod_mu;
%% This is only run for CFC11 for the 'unexpected emissions' scenario
if fugitive_emissions
    Fug_emiss = 13*10^3; %tonnes/yr since 2012 from Montzka et al. 2018
    ytmp1 = find(years == 2012);
    Prod_UL(ytmp1:end) = Fug_emiss/mean(DE_sim(end,:));
    AddtoUB = linspace(0,Prod_UL(ytmp1),13);
    % For 10 years prior to 2012, add linear increase of the parameter for
    % the upper limit
    Prod_UL(ytmp1-12:ytmp1-1) = Prod_UL(ytmp1-12:ytmp1-1)+AddtoUB(1:12);
    if meanSPARC
        OutputFileName = strcat(InputFileName,'_meanSparcLT_fugitive.mat');
    elseif Constant_LT
        OutputFileName = strcat(InputFileName,'_constantLT_fugitive.mat');
    end        
end
%% Observationally derived emissions
load('wmo2018.mat');
ytmp1 = find(wmo_yr == y1);
ytmp2 = find(wmo_yr == y2);
conc_obs = wmo_conc(ytmp1:ytmp2+1);

ObsDerivedEmiss = ppt_to_tonnes*(conc_obs(2:end)'-conc_obs(1:end-1)'.*exp(-1./LT(1:end-1)'));
ObsDerivedEmiss_max = ppt_to_tonnes*(conc_obs(2:end)'-conc_obs(1:end-1)'.*exp(-1./LT_min(1:end-1)'));
ObsDerivedEmiss_min = ppt_to_tonnes*(conc_obs(2:end)'-conc_obs(1:end-1)'.*exp(-1./LT_max(1:end-1)'));

% Saving mean of RF and DE and reported Prod values (named Prod_mu)
RF_mu = mean(RF_sim,2); 
DE_mu = mean(DE_sim,2); 

cd(HomeDir)
OutputFolder_FileName = strcat(InputFileName,'/Output/meanvars',OutputFileName);
save(OutputFolder_FileName,'ObsDerivedEmiss','Prod_mu','RF_mu','DE_mu','LT','conc_obs')


%% Run the Simulation Model
if fugitive_emissions
    [Bank,Emiss,Prod_prior,IND,rho1,rho3] = simulation_model(Prod_UL, Prod_LL, DE_sim, RF_sim, BankSize1,N,1);
else
    [Bank,Emiss,Prod_prior,IND,rho1] = simulation_model(Prod_UL, Prod_LL, DE_sim, RF_sim, BankSize1,N,0);
end

%% Calculate Top-down estimate of banks
Banks_TD(1,1) = BankSize1 + Prod_mu(1)- ObsDerivedEmiss(1); 
for t = 2:Nyears
    Banks_TD(t,1) = Banks_TD(t-1) + Prod_mu(t)- ObsDerivedEmiss(t);
end
    
% Emissions uncertainty samples, based on when observations are available.
ytmp1 = find(years == 1980);
ytmp2 = find(years == 2016);
emiss_diff = abs(mean(ObsDerivedEmiss_max(ytmp1:ytmp2)-ObsDerivedEmiss_min(ytmp1:ytmp2)));
emiss_sigma_prior = 2*max(emiss_diff,2*10^4)*repmat(betarnd(5,5,N,1),1,Nyears);

% Estimate rho_err, the autocorrelation of covariance matrix
% Ranges are based on previous tests
if fugitive_emissions
    rho_err = 0.7+0.3*betarnd(2,2,N,1);
elseif CFC113
    rho_err = betarnd(2,2,N,1);            
else
    rho_err = 0.5+0.5*betarnd(2,2,N,1);
end

ytmp1 = find(years == 1980);
ytmp2 = find(years == 2016);
N_tmp = length(ObsDerivedEmiss(ytmp1:end));
for ii = 1:N
    % Adding correlation to covariance matrix
    rho_tmp = rho_err(ii)*ones(N_tmp,N_tmp);
    exp_val = abs(repmat([1:N_tmp],N_tmp,1)-repmat([1:N_tmp]',1,N_tmp));
    Rhom = rho_tmp.^exp_val;
    Cov_matrix = Rhom.*(emiss_sigma_prior(ii,ytmp1:end)'*emiss_sigma_prior(ii,ytmp1:end));
    % mvnpdf is a matlab function that calculates the joint probability of Emiss given ObsDerivedEmiss
    likelihood_fcn(ii) = mvnpdf(Emiss(ii,ytmp1:end),ObsDerivedEmiss(ytmp1:end),Cov_matrix);
end

IR_vec = (1/sum(likelihood_fcn))*cumsum(likelihood_fcn); % Importance Ratio vector

figure; plot(IR_vec); % Check that resamples come from enough samples

% Resample from the priors based on the relative likelihood
indx = nan(Nresamps,1);
parfor ii = 1:Nresamps
    indx(ii) = find(IR_vec>rand,1);
end

ResampleIndex = indx;

RF_prior = RF_sim(:,IND); 
DE_prior = DE_sim(:,IND); 

OutputFolder_FileName = strcat(OutputFolderName,OutputFileName);     
save(OutputFolder_FileName,'RF_prior','DE_prior','Bank','Emiss','ResampleIndex','Prod_prior','Banks_TD','ObsDerivedEmiss','emiss_sigma_prior','rho_err','rho1','IND')