% Yong Wang -- 2015.07.17
% Thanks to the help of Salvalaglio Matteo
% 
clc
clear all
close all 

load('times.txt');

fileID = fopen('SD_Nsim.txt','w');
fprintf(fileID,'#  M  <mu>       SD_mu      <t_m>      SD_t_m    <P>    SD_P    <tau>      SD_tau \n')

nboot=20; 


for M=[10:2:length(times)]
disp(M)   
for i=1:nboot;
    disp(i)
    t=datasample(times,M,'Replace',true); %HERE M is the number of simulations you want to consider

    Hi = STP_noplot(t,min(t)./1E2,max(t).*1E2,1E5,0.05,10);
    %t=times(B(:,i));
    %Hi=STP_noplot(t,1E-6,max(t).*10,1E4,0.05)
    Pi(i)=Hi.pvalue_KS_statistic;
    TAUi(i)=Hi.tau;
    MUi(i)=Hi.mu;
    t_mi(i)=Hi.t_m;

end

ERROR_mu(M)=std(MUi);
MEAN_mu(M)=mean(MUi);
ERROR_t_m(M)=std(t_mi);
MEAN_t_m(M)=mean(t_mi);
ERROR_P(M)=std(Pi);
MEAN_P(M)=mean(Pi);
% that’s the standard deviation of the fitted time, its order of magnitude is typically consistent with the std(mu)
ERROR_Tau(M)=std(TAUi);
MEAN_Tau(M)=mean(TAUi);

formatSpec = '%4s is ( %10.5f +- %10.5f ) [time unit] \n';
fprintf('================================ \n');
fprintf('Number of simulations = %4d \n',M)
fprintf(formatSpec,'mu',MEAN_mu(M),ERROR_mu(M))
fprintf(formatSpec,'t_m',MEAN_t_m(M),ERROR_t_m(M))
fprintf(formatSpec,'P',MEAN_P(M),ERROR_P(M))
fprintf(formatSpec,'tau',MEAN_Tau(M),ERROR_Tau(M))
fprintf('================================ \n');

%fprintf(fileID,'================================ \n');
%fprintf(fileID,'Number of simulations = %4d \n',M)
%fprintf(fileID,formatSpec,'mu',MEAN_mu(M),ERROR_mu(M))
%fprintf(fileID,formatSpec,'t_m',MEAN_t_m(M),ERROR_t_m(M))
%fprintf(fileID,formatSpec,'P',MEAN_P(M),ERROR_P(M))
%fprintf(fileID,formatSpec,'tau',MEAN_Tau(M),ERROR_Tau(M))
%fprintf(fileID,'================================ \n');
fprintf(fileID,'%4d %10.5f %10.5f %10.5f %10.5f %6.4f %6.4f %10.5f %10.5f \n',M,MEAN_mu(M),ERROR_mu(M),MEAN_t_m(M),ERROR_t_m(M),MEAN_P(M),ERROR_P(M),MEAN_Tau(M),ERROR_Tau(M))

end

%H=STP_noplot(times,min(times)./10,max(times).*10,1E4,0.05)
H=STPtest(times,min(times)./10,max(times).*10,1E4,0.05)
fclose(fileID);
