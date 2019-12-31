%estimating Relese Fractions using AFEAS data and Ashford release fractions

clear all
HomeDir = '/Users/meganlickley/Dropbox (MIT)/Code';
str = strcat(HomeDir,'/CFC11/Input/AFEAS_cfc11production.mat');
load(str,'closecell','nonhermetic','open_aero','yr')
Yend = 2018;
Nyp = length(yr(end)+1:Yend);

N = 5*10^5;
Pu = 0.1*lognrnd(0,0.5,length(yr)+Nyp,N,3)+0.95;%Production Uncertainty

RFvals.aerosol(1,:) = betarnd(8,6,1,N);
RFvals.aerosol(2,:) = 1-RFvals.aerosol(1,:);

m = 0.0366; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFvals.closedcell = lognrnd(mu,sigma,1,N);

m = 0.07; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFvals.nonhermetic = lognrnd(mu,sigma,1,N);

m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFvals.nonhermeticLT = lognrnd(mu,sigma,1,N);
Bank(1,:) = zeros(1,N);

closecell1 = repmat(closecell(1),1,N).*Pu(1,:,1);
nonhermetic1 = repmat(nonhermetic(1),1,N).*Pu(1,:,2);
open_aero1 = repmat(open_aero(1),1,N).*Pu(1,:,3);

Prod_AFEAS(1,:) = closecell1+nonhermetic1+open_aero1;
DE(1,:) = RFvals.nonhermetic.*nonhermetic1+RFvals.aerosol(1,:).*open_aero1;
NH_bank(1,:) = (1-RFvals.nonhermetic).*nonhermetic1;
CC_bank(1,:) = closecell1;
Bank_In(1,:) = Prod_AFEAS(1,:)-DE(1,:);
for y = 2: length(yr)
    NH(y,:) = repmat(nonhermetic(y)-nonhermetic(y-1),1,N).*Pu(y,:,1);
    OA(y,:) = repmat(open_aero(y)-open_aero(y-1),1,N).*Pu(y,:,2);
    CC(y,:) = repmat(closecell(y)-closecell(y-1),1,N).*Pu(y,:,3);

    DE(y,:) = RFvals.nonhermetic.*NH(y,:)+RFvals.aerosol(1,:).*OA(y,:);
    Prod_AFEAS(y,:) = NH(y,:)+OA(y,:)+CC(y,:);

    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    NH_bank(y,:) = exp(-ones(1,N)./RFvals.nonhermeticLT).*NH_bank(y-1,:)+(1-RFvals.nonhermetic).*NH(y,:);
    CC_bank(y,:) = (1-RFvals.closedcell).*CC_bank(y-1,:)+CC(y,:);
    Bank_emiss(y,:) = RFvals.aerosol(2,:).*OA(y-1,:)+NH_bank(y-1,:).*(1-exp(-ones(1,N)./RFvals.nonhermeticLT))...
        + RFvals.closedcell.*CC_bank(y-1,:);
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

% For years after AFEAS production ends
for y = length(yr)+1:length(yr)+Nyp
    NH(y,:) = repmat(nonhermetic(end)-nonhermetic(end-1),1,N).*Pu(y,:,1);
    OA(y,:) = repmat(open_aero(end)-open_aero(end-1),1,N).*Pu(y,:,2);
    CC(y,:) = repmat(closecell(end)-closecell(end-1),1,N).*Pu(y,:,3);

    DE(y,:) = RFvals.nonhermetic.*NH(y,:)+RFvals.aerosol(1,:).*OA(y,:);
    Prod_AFEAS(y,:) = NH(y,:)+OA(y,:)+CC(y,:);

    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    NH_bank(y,:) = exp(-ones(1,N)./RFvals.nonhermeticLT).*NH_bank(y-1,:)+(1-RFvals.nonhermetic).*NH(y,:);
    CC_bank(y,:) = (1-RFvals.closedcell).*CC_bank(y-1,:)+CC(y,:);
    Bank_emiss(y,:) = RFvals.aerosol(2,:).*OA(y-1,:)+NH_bank(y-1,:).*(1-exp(-ones(1,N)./RFvals.nonhermeticLT))...
        + RFvals.closedcell.*CC_bank(y-1,:);
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

RF = Bank_emiss(2:end,:)./Bank(1:end-1,:);
RF = [NaN(1,N); RF];
DE = DE./Prod_AFEAS;

yr = yr(1):Yend;
str = strcat(HomeDir,'/CFC11/Input/InferredRFandDE.mat');

save(str,'RF','yr', 'DE','NH_bank','CC_bank')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%  Fugitive emission scenario %%%%%%%%%%%%%%%%%%%

clearvars -except HomeDir  N Yend

str = strcat(HomeDir,'/CFC11/Input/AFEAS_cfc11production.mat');
load(str,'closecell','nonhermetic','open_aero','yr')

Nyp = length(yr(end)+1:Yend);

Pu = 0.1*lognrnd(0,0.5,length(yr)+Nyp,N,3)+0.95;%Production Uncertainty
RFvals.aerosol(1,:) = betarnd(8,6,1,N);
RFvals.aerosol(2,:) = 1-RFvals.aerosol(1,:);

m = 0.0366; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFvals.closedcell = lognrnd(mu,sigma,1,N);

m = 0.07; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFvals.nonhermetic = lognrnd(mu,sigma,1,N);

m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RFvals.nonhermeticLT = lognrnd(mu,sigma,1,N);


closecell1 = repmat(closecell(1),1,N).*Pu(1,:,1);
nonhermetic1 = repmat(nonhermetic(1),1,N).*Pu(1,:,2);
open_aero1 = repmat(open_aero(1),1,N).*Pu(1,:,3);
Bank(1,:) = zeros(1,N);

Prod_AFEAS(1,:) = closecell1+nonhermetic1+open_aero1;
DE(1,:) = RFvals.nonhermetic.*nonhermetic1+RFvals.aerosol(1,:).*open_aero1;
NH_bank(1,:) = (1-RFvals.nonhermetic).*nonhermetic1;
CC_bank(1,:) = closecell1;
Bank_In(1,:) = Prod_AFEAS(1,:)-DE(1,:);
for y = 2: length(yr)
    NH(y,:) = repmat(nonhermetic(y)-nonhermetic(y-1),1,N).*Pu(y,:,1);
    OA(y,:) = repmat(open_aero(y)-open_aero(y-1),1,N).*Pu(y,:,2);
    CC(y,:) = repmat(closecell(y)-closecell(y-1),1,N).*Pu(y,:,3);

    DE(y,:) = RFvals.nonhermetic.*NH(y,:)+RFvals.aerosol(1,:).*OA(y,:);
    Prod_AFEAS(y,:) = NH(y,:)+OA(y,:)+CC(y,:);

    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    NH_bank(y,:) = exp(-ones(1,N)./RFvals.nonhermeticLT).*NH_bank(y-1,:)+(1-RFvals.nonhermetic).*NH(y,:);
    CC_bank(y,:) = (1-RFvals.closedcell).*CC_bank(y-1,:)+CC(y,:);
    Bank_emiss(y,:) = RFvals.aerosol(2,:).*OA(y-1,:)+NH_bank(y-1,:).*(1-exp(-ones(1,N)./RFvals.nonhermeticLT))...
        + RFvals.closedcell.*CC_bank(y-1,:);
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

tmp = rand(N,3);

NHextra = tmp(:,1)./sum(tmp,2);
OAextra = tmp(:,2)./sum(tmp,2);
CCextra = tmp(:,3)./sum(tmp,2);

Prod_2012 = 71; % From Montzka et al. 2018 estimate
FugitiveADD = linspace(0,Prod_2012,13);
FugitiveADD(end+1:18) = Prod_2012;
Prod_randvals = rand(N,1);

for y = length(yr)+1:length(yr)+Nyp
    NH(y,:) = repmat(nonhermetic(end)-nonhermetic(end-1),1,N).*Pu(y,:,1)+FugitiveADD(y-length(yr))*NHextra'.*Prod_randvals';
    OA(y,:) = repmat(open_aero(end)-open_aero(end-1),1,N).*Pu(y,:,2)+FugitiveADD(y-length(yr))*OAextra'.*Prod_randvals';
    CC(y,:) = repmat(closecell(end)-closecell(end-1),1,N).*Pu(y,:,3)+FugitiveADD(y-length(yr))*CCextra'.*Prod_randvals';

    DE(y,:) = RFvals.nonhermetic.*NH(y,:)+RFvals.aerosol(1,:).*OA(y,:);
    Prod_AFEAS(y,:) = NH(y,:)+OA(y,:)+CC(y,:);

    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    NH_bank(y,:) = exp(-ones(1,N)./RFvals.nonhermeticLT).*NH_bank(y-1,:)+(1-RFvals.nonhermetic).*NH(y,:);
    CC_bank(y,:) = (1-RFvals.closedcell).*CC_bank(y-1,:)+CC(y,:);
    Bank_emiss(y,:) = RFvals.aerosol(2,:).*OA(y-1,:)+NH_bank(y-1,:).*(1-exp(-ones(1,N)./RFvals.nonhermeticLT))...
        + RFvals.closedcell.*CC_bank(y-1,:);
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

RF = Bank_emiss(2:end,:)./Bank(1:end-1,:);
RF = [NaN(1,N); RF];
DE = DE./Prod_AFEAS;

yr = yr(1):Yend;
str = strcat(HomeDir,'/CFC11/Input/Fugitive_InferredRFandDE.mat');

save(str,'RF','yr', 'DE')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CFC-12 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars -except HomeDir  N Yend

str = strcat(HomeDir,'/CFC12/Input/AFEAS_cfc12production.mat');
load(str,'closecell','nonhermetic','open_aero','refrigeration','yr')

Nyp = length(yr(end)+1:Yend);

Pu = 0.1*lognrnd(0,0.5,length(yr)+Nyp,N,4)+0.95;

% model aerosol and opencell foam together bc of AFEAS grouping
RF12.aerosol(1,:) = betarnd(8,6,1,N);
RF12.aerosol(2,:) = 1-RF12.aerosol(1,:);

RF12.closedcell(1,:) = betarnd(12,12,1,N);
RF12.closedcell(2,:) = 1-RF12.closedcell(1,:);

m = 0.07; v = (0.5*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF12.nonhermetic = lognrnd(mu,sigma,1,N);%0.07*DEu(5,:);

m = 10; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF12.nonhermeticLT = lognrnd(mu,sigma,1,N);

m = 0.02; v = (0.7*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF12.hermetic = lognrnd(mu,sigma,1,N); %0.02*DEu(7,:);

m = 20; v = (0.20*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF12.hermeticLT = lognrnd(mu,sigma,1,N);

closecell1 = repmat(closecell(1),1,N).*Pu(1,:,1);
nonhermetic1 = repmat(nonhermetic(1),1,N).*Pu(1,:,2);
open_aero1 = repmat(open_aero(1),1,N).*Pu(1,:,3);
refrigeration1 =  repmat(refrigeration(1),1,N).*Pu(1,:,4);

Prod_AFEAS(1,:) = closecell1+nonhermetic1+open_aero1+refrigeration1;
DE(1,:) = RF12.nonhermetic.*nonhermetic1+RF12.aerosol(1,:).*open_aero1+RF12.hermetic.*refrigeration1+RF12.closedcell(1,:).*closecell1;
NH_bank(1,:) = (ones(1,N)-RF12.nonhermetic).*nonhermetic1;
HR_bank(1,:) = (ones(1,N)-RF12.hermetic).*refrigeration1;
Bank_In(1,:) = Prod_AFEAS(1,:)-DE(1,:);
Bank(1,:) = Bank_In(1,:);
for y = 2: length(yr)
    NH(y,:) = repmat((nonhermetic(y)-nonhermetic(y-1)),1,N).*Pu(y,:,1);
    HR(y,:) = repmat((refrigeration(y) - refrigeration(y-1)),1,N).*Pu(y,:,2);
    OA(y,:) = repmat((open_aero(y)-open_aero(y-1)),1,N).*Pu(y,:,3);
    CC(y,:) = repmat((closecell(y)-closecell(y-1)),1,N).*Pu(y,:,4);
    
    DE(y,:) = RF12.nonhermetic.*NH(y,:)+RF12.aerosol(1,:).*OA(y,:)+RF12.hermetic.*HR(y,:)...
        +RF12.closedcell(1,:).*CC(y,:);
    Prod_AFEAS(y,:) = NH(y,:)+OA(y,:)+CC(y,:)+HR(y,:);
    
    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    NH_bank(y,:) = exp(-ones(1,N)./RF12.nonhermeticLT).*NH_bank(y-1,:)+(1-RF12.nonhermetic).*NH(y,:);
    HR_bank(y,:) = exp(-ones(1,N)./RF12.hermeticLT).*HR_bank(y-1,:)+(1-RF12.hermetic).*HR(y,:);
    
    Bank_emiss(y,:) = RF12.aerosol(2,:).*OA(y-1,:)+NH_bank(y-1,:).*(1-exp(-ones(1,N)./RF12.nonhermeticLT))...
        + RF12.closedcell(2,:).*CC(y-1,:)+HR_bank(y-1,:).*(1-exp(-ones(1,N)./RF12.hermeticLT));
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

for y = length(yr)+1: length(yr)+Nyp
    NH(y,:) = repmat((nonhermetic(end)-nonhermetic(end-1)),1,N).*Pu(y,:,1);
    HR(y,:) = repmat((refrigeration(end) - refrigeration(end-1)),1,N).*Pu(y,:,2);
    OA(y,:) = repmat((open_aero(end)-open_aero(end-1)),1,N).*Pu(y,:,3);
    CC(y,:) = repmat((closecell(end)-closecell(end-1)),1,N).*Pu(y,:,4);
    
    DE(y,:) = RF12.nonhermetic.*NH(y,:)+RF12.aerosol(1,:).*OA(y,:)+RF12.hermetic.*HR(y,:)...
        +RF12.closedcell(1,:).*CC(y,:);
    Prod_AFEAS(y,:) = NH(y,:)+OA(y,:)+CC(y,:)+HR(y,:);
    
    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:); 
    NH_bank(y,:) = exp(-ones(1,N)./RF12.nonhermeticLT).*NH_bank(y-1,:)+(1-RF12.nonhermetic).*NH(y,:);
    HR_bank(y,:) = exp(-ones(1,N)./RF12.hermeticLT).*HR_bank(y-1,:)+(1-RF12.hermetic).*HR(y,:);
    
    Bank_emiss(y,:) = RF12.aerosol(2,:).*OA(y-1,:)+NH_bank(y-1,:).*(1-exp(-ones(1,N)./RF12.nonhermeticLT))...
        + RF12.closedcell(2,:).*CC(y-1,:)+HR_bank(y-1,:).*(1-exp(-ones(1,N)./RF12.hermeticLT));
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

RF = Bank_emiss(2:end,:)./Bank(1:end-1,:);
RF = [NaN(1,N);RF];

DE = DE./Prod_AFEAS;

yr = yr(1):Yend;
str = strcat(HomeDir,'/CFC12/Input/InferredRFandDE.mat');
save(str,'RF','yr', 'DE','HR_bank','NH_bank')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CFC-113 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars -except HomeDir N Yend

str = strcat(HomeDir,'/CFC113/Input/AFEAS_cfc113production.mat');
load(str,'longbank','shortbank','yr')

Nyp = length(yr(end)+1:Yend);

Pu = 0.1*lognrnd(0,0.5,length(yr)+Nyp,N,4)+0.95;

RF113.SB(1,:) = betarnd(12,12,1,N);
RF113.SB(2,:) = 1-RF113.SB(1,:);

m = 0.02; v = (0.50*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF113.LB = lognrnd(mu,sigma,1,N); %0.02*DEu(2,:);

m = 20; v = (0.2*m)^2;
mu = log(m^2/sqrt(v+m^2));
sigma = sqrt(log(v/m^2+1));
RF113.LB_LT = lognrnd(mu,sigma,1,N); 

shortbank1 = repmat(shortbank(1),1,N).*Pu(1,:,1);
longbank1 = repmat(longbank(1),1,N).*Pu(1,:,2);
Prod_AFEAS(1,:) = shortbank1+longbank1;
DE(1,:) = RF113.SB(1,:).*shortbank1+RF113.LB.*longbank1;
LB_bank(1,:) = (1-RF113.LB).*longbank1;
Bank_In(1,:) = Prod_AFEAS(1,:)-DE(1,:);
Bank(1,:) = Bank_In(1,:);

for y = 2: length(yr)
    SB(y,:) = repmat(shortbank(y) - shortbank(y-1),1,N).*Pu(y,:,1);
    LB(y,:) = repmat(longbank(y) - longbank(y-1),1,N).*Pu(y,:,2);
    
    DE(y,:) = RF113.SB(1,:).*SB(y,:)+RF113.LB.*LB(y,:);
    Prod_AFEAS(y,:) = SB(y,:)+LB(y,:);
    
    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    
    LB_bank(y,:) = exp(-ones(1,N)./RF113.LB_LT).*LB_bank(y-1,:)+(1-RF113.LB).*LB(y,:);
    
    Bank_emiss(y,:) = RF113.SB(2,:).*SB(y-1,:)+LB_bank(y-1,:).*(1-exp(-ones(1,N)./RF113.LB_LT));
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

for y = length(yr)+1: length(yr)+Nyp
    SB(y,:) = repmat(shortbank(end) - shortbank(end-1),1,N).*Pu(y,:,1);
    LB(y,:) = repmat(longbank(end) - longbank(end-1),1,N).*Pu(y,:,2);
    
    DE(y,:) = RF113.SB(1,:).*SB(y,:)+RF113.LB.*LB(y,:);
    Prod_AFEAS(y,:) = SB(y,:)+LB(y,:);
    
    Bank_In(y,:) = Prod_AFEAS(y,:)-DE(y,:);
    
    LB_bank(y,:) = exp(-ones(1,N)./RF113.LB_LT).*LB_bank(y-1,:)+(1-RF113.LB).*LB(y,:);
    
    Bank_emiss(y,:) = RF113.SB(2,:).*SB(y-1,:)+LB_bank(y-1,:).*(1-exp(-ones(1,N)./RF113.LB_LT));
    Bank(y,:) = Bank(y-1,:) + Bank_In(y,:) - Bank_emiss(y,:);
end

RF = Bank_emiss(2:end,:)./Bank(1:end-1,:);
RF = [NaN(1,N);RF];

DE = DE./Prod_AFEAS;

yr = yr(1):Yend;
str = strcat(HomeDir,'/CFC113/Input/InferredRFandDE.mat');
save(str,'RF','yr', 'DE')
