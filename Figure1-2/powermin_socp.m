function [feasible,Wsolution] = powermin_socp(H, L, N1, A_set, P , r, delta)
%The implementation utilizes and requires CVX: http://cvxr.com/

%INPUT:
%H          = Kt*Nt x Kr matrix with row index for transmit and column
%             index receiver antennas
%L         =# RRHs
%N1       =# antennas of each RRH
%A_set   = RRHs set
%P          =1 x L vector Limits of the L power constraints
%r          = 1 x K vector with SINR constraints for all users.

%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = N x k matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.

K = size(H,2); %Number of users

cvx_begin quiet
      variable W(L*N1,K) complex;   %Variable for N x K beamforming matrix
      variable t(length(A_set),1); %Auxiliary variable for power constraints
minimize (sum(t))
     subject to
     for k=1:K        %%%%%%%%%QoS constraints
         norm([H(:,k)'*W, delta])<=sqrt(1+1/r)*real(H(:,k)'*W(:,k));
     end
     for l=1:length(A_set)    %%%Active RRHs: Transmit Power Constraints
         norm(W(N1*(A_set(l)-1)+1:N1*A_set(l),:),'fro')<=t(l);
         t(l)<=sqrt(P(A_set(l)));
         t(l)>=0; 
     end
 cvx_end
     
     %Analyze result and prepare the output variables.
     if  strfind(cvx_status,'Solved') 
         feasible=true;
         Wsolution=W;
     else
         feasible=false;
         Wsolution=[];
     end


