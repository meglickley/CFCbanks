function [Bank,Emiss,Prod,indx,rho1,rho3] = forward_propagation(Prod_UL, Prod_LL, DE_sim, RF_sim, BankSize1,N,fugitive)

t = 1;
Nyears = length(Prod_LL);

Nyears1 = 34; % Change point between AFEAS and UNEP data

rho1(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
rho_tmp = repmat(rho1,Nyears1,Nyears1,1);
exp_val = repmat(abs(repmat([1:Nyears1],Nyears1,1)-repmat([1:Nyears1]',1,Nyears1)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyears1), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
Prod(:,1:Nyears1) = (repmat(0.2*Prod_LL(1:Nyears1),N,1)).*logninv(Um,zeros(N,Nyears1),0.5*ones(N,Nyears1))+0.95*repmat(Prod_LL(1:Nyears1),N,1);
clear Um

Nyears2 = Nyears - Nyears1; 
rho2 = rho1;
rho_tmp = repmat(rho2,Nyears2,Nyears2,1);

exp_val = repmat(abs(repmat([1:Nyears2],Nyears2,1)-repmat([1:Nyears2]',1,Nyears2)),1,1,N);
Rhom = rho_tmp.^exp_val;
clear rho_tmp exp_val
Zm = mvnrnd(zeros(N,Nyears2), Rhom, N);
clear Rhom
Um = normcdf(Zm,0,1);
clear Zm
Prod(:,Nyears1+1:Nyears) = (repmat(0.1*Prod_LL(Nyears1+1:Nyears),N,1)).*logninv(Um,zeros(N,Nyears2),0.5*ones(N,Nyears2))+0.95*repmat(Prod_LL(Nyears1+1:Nyears),N,1);
clear Um

if fugitive
    Nyears3 = Nyears-46;
    rho3(1,1,:) = 0.5+0.5*betarnd(2,2,N,1);
    rho_tmp = repmat(rho3,Nyears3,Nyears3,1);

    exp_val = repmat(abs(repmat([1:Nyears3],Nyears3,1)-repmat([1:Nyears3]',1,Nyears3)),1,1,N);
    Rhom = rho_tmp.^exp_val;
    clear rho_tmp exp_val
    Zm = mvnrnd(zeros(N,Nyears3), Rhom, N);
    clear Rhom
    Um = normcdf(Zm,0,1);
    clear Zm
    Prod(:,46+1:Nyears) = (repmat(Prod_UL(46+1:Nyears)-Prod_LL(46+1:Nyears),N,1)).*betainv(Um,ones(N,Nyears3),ones(N,Nyears3))+repmat(Prod_LL(46+1:Nyears),N,1);
    clear Um
end

%%

indx = randi(size(DE_sim,2),N,1);
DEsamps = DE_sim(t,indx)';
RFsamps = RF_sim(t,indx)';
Prodsamps = Prod(:,t);
[bank_samps,emiss_samps] = Bank_Emiss(Prodsamps,RFsamps,DEsamps,BankSize1*ones(size(Prodsamps)));
Bank(:,t) = bank_samps;
Emiss(:,t) = emiss_samps;
for t = 2:Nyears
    DEsamps = DE_sim(t,indx)';
    RFsamps = RF_sim(t,indx)';
    Prodsamps = Prod(:,t);
    [bank_samps,emiss_samps] = Bank_Emiss(Prodsamps,RFsamps,DEsamps,bank_samps);
    Bank(:,t) = bank_samps;
    Emiss(:,t) = emiss_samps; 
end