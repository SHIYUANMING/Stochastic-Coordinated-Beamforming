%%%%%%%%%%Main Setup%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all; warning off all;
cvx_solver mosek;
L=5; K=3; N1=1; %%% #RRH:L, #user: K, #antenna: N1 
delta=1;   %normized noise variance
P=10^(0)*ones(1,L); % P: power constraints
Q=3;   % QoS requirements  
r=(10^(Q/10));%%%%all the MU has the same QoS requirements 
epsilon=0.1; %%%Outage Probability
RRH_set=[1:L];   %A_set: active RRH set, 
R=3; % # Channel Coefficients that Can not be Obtained for Each MU
TT=100;  % # Iterations for the Algorithms
tauu=0.01;  % estimation errors
S1=308;

  load('D.mat');
  load('H.mat'); 

[Omega1, Omega2]=CompressiveCSI(D, R);   %%%%% Determine the Set Omega;
% 
% for tt=1:TT  %%Generate samples for Scenario Approach
% H_samples_Scenario_temp=samples(H, D, Omega1, Omega2, N1, S1, tauu); %Generate S1 Samples for S
% H_samples_Scenario(:,:,:,tt)=H_samples_Scenario_temp;
% end

load('H_Scenario.mat');

for tt=1:TT   %%%Iteration Numbers
H_samples_Scenario=H_Scenario(:,:,:,tt); 
%%%%%%%%%%%%%%%Scenario Approach%%%%%%%%%%%
[feasible_Scenario,W_Scenario]=powermin_socp_iterative(H_samples_Scenario, L, N1, RRH_set, P, r, delta, S1);
 if  feasible_Scenario==1    %feasilbe
     P_Scenario=norm(W_Scenario,'fro')^2;  %%total transmit power
 else
     P_Scenario=10^20;
 end
TotalPower_Scenario(tt)=P_Scenario;
W_Sceanrio_recoder(:,:,tt)=W_Scenario;
end

plot([0:TT-1],10*log10(TotalPower_Scenario.*1000),'b-o','LineWidth',2.5, 'MarkerSize',10); %Bernstein Approximation
hold on;

h=legend('Benchmark: Full CSI', 'Stochastic DC Programming','Fixed One','Fixed Two', 'fontsize',12,'fontweight','b','fontname','helvetica');
xlabel('Iteration','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Total Transmit Power [dBm]','fontsize',14,'fontweight','b','fontname','helvetica');