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
TT=3;  % # Iterations for the Algorithms
tauu=0.01;  % estimation errors
S_DC=1000;

  load('D.mat');
  load('H.mat'); 

[Omega1, Omega2]=CompressiveCSI(D, R);   %%%%% Determine the Set Omega; 
% for tt=1:TT  %%Generate samples for Scenario Approach
% H_samples_DC_temp=samples(H, D, Omega1, Omega2, N1, S_DC, tauu); %Generate S1 Samples for S
% H_DC(:,:,:,tt)=H_samples_DC_temp;
% end
% 
% save('H_DC.mat');

load('H_DC.mat');
load('W_Bernstein.mat');

for tt=1:TT   %%%Iteration Numbers
H_samples_DC=H_DC(:,:,:,tt); 

%%%%%%%%%%%%Stochastic DC Programming%%%%%%%%%%%
%%
%%%%%%%%Initial Solution%%%%%%%%%%%%%%%%%%%%%%%%%
feasible_Initial=1; %Initial Solution
W_iteration=W_Bernstein;  %recorder j-th iteration

P_iteration_1=10^10;  
P_iteration_2=0;

if feasible_Initial==1  %stop critria
iteration_eps=10^20;
kappa=10^20;
else
    iteration_eps=0;
end
DC_iteration=0;
while iteration_eps>10^(-4)
  [feasible_iteration,W_temp, kappa_temp, cvx_opt]=powermin_DC_iteration(H_samples_DC, W_iteration, S_DC, K, L, N1, r, delta, RRH_set, P, epsilon);
  DC_iteration=DC_iteration+1;

 if  feasible_iteration==1    %feasilbe
     kappa=kappa_temp;
     
     W_iteration=W_temp;
     P_iteration_2=cvx_opt;
     P_DC_2=norm(W_iteration,'fro')^2;  %%total transmit power
     W_DC_2=W_iteration;
     
     iteration_eps=abs(P_iteration_1-P_iteration_2);
     
     P_iteration_1=P_iteration_2;
     P_DC_1=P_DC_2;
     W_DC_1=W_DC_2;
     
 else %infeasible
   P_iteration_2=10^20;
   iteration_eps=0;
   kappa=0;
 end
 
end

TotalPower_DC(tt)=P_DC_1;
W_DC_recoder(:,:,tt)=W_DC_1;
kappa_recoder(tt)=kappa;

end


plot([0:length(TotalPower_DC)-1],10*log10(TotalPower_DC.*1000),'b-o','LineWidth',2.5, 'MarkerSize',10); %Bernstein Approximation
hold on;

h=legend('DC', 'fontsize',12,'fontweight','b','fontname','helvetica');
xlabel('Iteration','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Total Transmit Power [dBm]','fontsize',14,'fontweight','b','fontname','helvetica');

MySendMail;