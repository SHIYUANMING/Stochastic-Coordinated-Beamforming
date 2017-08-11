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
TT=41;  % # Iterations for the Algorithms
tauu=0.01;  %estimation errors

kappa_index=0;
W_cvx_temp=0;

S2=1000; % # Samples for the Stochastic DC Programming

%%%%%%%%%%%%%%%%Channel Realization%%%%%%%%%%%%%%%%%%%%%%%%
%  U_position=800.*(rand(2,K)-0.5);  %% user positions
%  B_position=800.*(rand(2,L)-0.5);   %%RRH positions
% for k=1:K
%     for l=1:L
%                  d=(norm(B_position(:,l)-U_position(:,k))+10);
%                  D(l,k)=4.4*10^(5)/(d^(1.88)*10^(normrnd(0,6.3)/20));
%                  H(N1*(l-1)+1:N1*l,k)=D(l,k)*(normrnd(0,1/sqrt(2),N1,1)+i*normrnd(0,1/sqrt(2),N1,1));  %%%nosie normalized to 1
%     end
% end
  load('D.mat');
  load('H.mat'); 

[Omega1, Omega2]=CompressiveCSI(D, R);   %%%%% Determine the Set Omega;
%H_samples_DC=samples(H, D, Omega1, Omega2, N1, S2, tauu); %Generate S2 Samples for Stochastic DC Programming
load('H_samples_DC.mat');

for tt=1:TT   %%%Iteration Numbers
%%%%%%%%%%%%%%%Bernstein Approximation%%%%%%%%%%%
if tt==1  %%%%%Only Solve Once
%[feasible_Bernstein,W_Bernstein] = powermin_Bernstein(H, D, Omega1, Omega2, tauu, L, K, N1, RRH_set, P, r, delta, epsilon/K);
%powermin_Bernstein_rankone(W_Bernstein, H, D, Omega1, Omega2, tauu, L, K, N1, RRH_set, P, r, delta, epsilon/K)

  load('W_Bernstein.mat');
  feasible_Bernstein=1;

if  feasible_Bernstein==1    %feasilbe
  P_Bernstein=norm(W_Bernstein,'fro')^2;  %%total transmit power
 else
   P_Bernstein=10^20;
 end
end
TotalPower_Bernstein(tt)=P_Bernstein;

%%%%%%%%%%%%%%%%Stochastic DC Programming%%%%%%%%%%%%%%
if tt==1  %initial solution
 W_iteration=W_Bernstein;
 feasible_iteration=feasible_Bernstein;
else
  [feasible_iteration,W_temp, kappa, W_cvx_temp]=powermin_DC_iteration(H_samples_DC, W_iteration, S2, K, L, N1, r, delta, RRH_set, P, epsilon);
  W_iteration=W_temp;
  kappa_index(tt)=kappa;
  W_cvx(tt)=W_cvx_temp;
  if tt>=2
      W_cvx(tt)-W_cvx(tt-1)
  end
end
 if  feasible_iteration==1    %feasilbe
  P_DC=norm(W_iteration,'fro')^2;  %%total transmit power
 end
TotalPower_DC(tt)=P_DC;
W_DC(:,:,tt)=W_iteration;

end

plot([0:length(TotalPower_Bernstein)-1],10*log10(TotalPower_Bernstein.*1000),'b-o','LineWidth',2.5, 'MarkerSize',10); %Bernstein Approximation
hold on;
plot([0:length(TotalPower_DC)-1],10*log10(TotalPower_DC.*1000),'g--', 'LineWidth',2.5, 'MarkerSize',8); %Stochastic DC Programming: Bernstein Initial
hold on;

h=legend('Benchmark: Full CSI', 'Stochastic DC Programming','Fixed Two', 'fontsize',12,'fontweight','b','fontname','helvetica');
xlabel('Iteration','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Total Transmit Power [dBm]','fontsize',14,'fontweight','b','fontname','helvetica');

MySendMail;