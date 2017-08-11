function [feasible,Wsolution,kappa_record, W_cvx] = powermin_DC_iteration(H_samples, W_iteration, M, K, L, N1, r, delta, A_set, P, epsilon)
%The implementation utilizes and requires CVX: http://cvxr.com/

%INPUT:
%H_samples       = Kt*Nt x Kr matrix with row index for transmit and column
%                          index receiver antennas, M# samples
%W_iteration       =# Beamforming matrix at the particular iteration
%M                    =# Samples 
%K                     =# Users
%L                     =# RRHs
%N1                   =# Antennas
%Mu                  = Parameter for conservative approximation
%r                     = QoS requirements
%delta               =Noise covariance 
%A_set              = RRHs set
%P                    =1 x L vector Limits of the L power constraints
%epsilon           =Outage constraint 


%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = N x k matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.

U_estimator=U_sampling(H_samples, W_iteration, M, K, 0, delta, r); %%%% 1x 1 scalar
Gradient_estimator=Gradient_sampling(H_samples, W_iteration, M, K, L, N1, delta, r); % L*N1 x 1 matrix
H=H_samples;

cvx_begin quiet
      variable W(L*N1,K) complex;   %Variable for N x K beamforming matrix
      variable y(M,1);  %Slack variables
      variable kappa
minimize norm(W,'fro')
     subject to
     kappa>0;
%%%%%%%%%linearization Constraint%%%%%%%%%%%%%%%
 (1/M)*sum(y)-real(U_estimator)-2*real((vec(W-W_iteration))'*Gradient_estimator)<=kappa*epsilon;   
 
 %%%%%%%%%%%%%Constraint (2) %%%%%%%%%%%%%%%%
for m=1:M
    for k=1:K  %%%Calculate S_k, k=1,...K
        temp=diag(H(:,:,m)'*W);
        temp1=H(:,:,m)'*W;
        if k==1
         real(temp1(k,[2:K])*temp1(k,[2:K])')+(1/r)*real(temp([2:K])'*temp([2:K]))+delta^2+kappa<=y(m);
        else
        real(temp1(k,[1:k-1,k+1:K])*temp1(k,[1:k-1,k+1:K])')+(1/r)*real(temp([1:k-1,k+1:K])'*temp([1:k-1,k+1:K]))+delta^2+kappa<=y(m);
        end
    end
    (1/r)*real(temp'*temp)<=y(m);
    y(m)>=0;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     for l=1:length(A_set)    %%%Active RRHs: Transmit Power Constraints
         norm(W(N1*(A_set(l)-1)+1:N1*A_set(l),:),'fro')<=sqrt(P(A_set(l)));
     end
 cvx_end
     
     %Analyze result and prepare the output variables.
     if  strfind(cvx_status,'Solved') 
         feasible=true;
         Wsolution=W;
         kappa_record=kappa;
         W_cvx=cvx_optval;
     else
         feasible=false;
         Wsolution=[];
         kappa_record=[];
         W_cvx=[];
     end