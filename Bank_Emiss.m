function [b,e] = Bank_Emiss_Conc(ProdIn,RFIn,DEIn,BanksIn,LT)
if size(DEIn,1) == 1
    b = ProdIn-DEIn*ProdIn+(ones(size(RFIn))-RFIn).*BanksIn;
    e = RFIn.*BanksIn+DEIn*ProdIn;
else
    b = ProdIn-DEIn.*ProdIn+(ones(size(RFIn))-RFIn).*BanksIn;
    e = RFIn.*BanksIn+DEIn.*ProdIn;
end
end