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
tauu=0.01;  % estimation errors

kappa_index=0;
W_cvx=0;
W_fixed=0;

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
H_samples_DC=samples(H, D, Omega1, Omega2, N1, S2, tauu); %Generate S2 Samples for Stochastic DC Programming

%load('H_samples_DC.mat');

for tt=1:TT   %%%Iteration Numbers
%%%%%%%%%%%%%%%Bernstein Approximation%%%%%%%%%%%
if tt==1  %%%%%Only Solve Once
[feasible_Bernstein,W_Bernstein] = powermin_Bernstein(H, D, Omega1, Omega2, tauu, L, K, N1, RRH_set, P, r, delta, epsilon/K);
powermin_Bernstein_rankone(W_Bernstein, H, D, Omega1, Omega2, tauu, L, K, N1, RRH_set, P, r, delta, epsilon/K)

%   load('W_Bernstein.mat');
%   feasible_Bernstein=1;

if  feasible_Bernstein==1    %feasilbe
  P_Bernstein=norm(W_Bernstein,'fro')^2;  %%total transmit power
 else
   P_Bernstein=10^20;
 end
end
TotalPower_Bernstein(tt)=P_Bernstein;


% %%%%%%%%%%%%%%%%Stochastic DC Programming%%%%%%%%%%%%%%%
if tt==1  %initial solution
 W_fixed=W_Bernstein;
 feasible_fixed=feasible_Bernstein;
else
  [feasible_fixed,W_temp, W_cvx_fixed_temp]=powermin_DC_iteration_fixed(0.01,H_samples_DC, W_fixed, S2, K, L, N1, r, delta, RRH_set, P, epsilon);
  W_fixed=W_temp;
  W_cvx_fixed(tt)=W_cvx_fixed_temp;
  if tt>=2
      W_cvx_fixed(tt-1)-W_cvx_fixed(tt)
  end
end
 if  feasible_fixed==1    %feasilbe
  P_fixed=norm(W_fixed,'fro')^2;  %%total transmit power 
 end
TotalPower_DC_fixed(tt)=P_fixed;
W_DC_fixed(:,:,tt)=W_fixed;
end

plot([0:TT-1],10*log10(TotalPower_Bernstein.*1000),'b-o','LineWidth',2.5, 'MarkerSize',10); %Bernstein Approximation
hold on;
plot([0:TT-1],10*log10(TotalPower_DC_fixed.*1000),'--','Color', [0,0,128]./256, 'LineWidth',2.5, 'MarkerSize',8); %Stochastic DC Programming: Bernstein Initial
hold on;

h=legend('Benchmark: Full CSI', 'Stochastic DC Programming','Fixed One','Fixed Two', 'fontsize',12,'fontweight','b','fontname','helvetica');
xlabel('Iteration','fontsize',14,'fontweight','b','fontname','helvetica');
ylabel('Total Transmit Power [dBm]','fontsize',14,'fontweight','b','fontname','helvetica');

MySendMail;