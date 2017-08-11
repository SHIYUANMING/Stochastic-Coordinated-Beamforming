function [feasible,Wsolution] = powermin_Bernstein_rankone(Wsolution_temp, H, D, Omega1, Omega2, tauu, L, K, N1, A_set, P, r, delta, epsilon)
%The implementation utilizes and requires CVX: http://cvxr.com/
%
%INPUT:
%H          = Kt*Nt x Kr matrix with row index for transmit and column
%D          =Average Network Channel Gain (Naive Scheme or Certainty Equivalent Scheme)
%Omega:  |Omega| x 2 matrix Positions of CSI that can only obtain statistical CSI

%L         =# RRHs
%K        =# MUs
%N1       =# antennas of each RRH
%A_set   = active RRHs set

%P          =1 x L vector Limits of the L power constraints
%r          = 1 x K vector with SINR constraints for all users.

%epsilon  = The outage constraints

%OUTPUT:
%feasible  = This variable is feasible=true if the feasibility problem is
%            feasible. Otherwise we have feasible=false.
%Wsolution = N x k matrix with beamforming achieved by the power
%            minimization problem. This matrix is empty if this problem is
%            infeasible.

%%%%%%%%%Caculate the parameters omega=CC for instantaneous CSI;  Lambda=DD for statistical CSI%%%%%%%%%%%
NN=size(Omega1,1); 
CC=H;   DD=zeros(L*N1, K); %Initialize 

for n=1:NN
CC(Omega1(n,1),Omega1(n,2))=0;   %%%%%set the unobserved CSI to zero
DD(Omega1(n,1),Omega1(n,2))=D(Omega1(n,1), Omega1(n,2));
end

NN2=size(Omega2,1);
for n=1:NN2
CC(Omega2(n,1),Omega2(n,2))=sqrt(1-tauu^2)*H(Omega2(n,1),Omega2(n,2));   %%%%%set the unobserved CSI to zero
DD(Omega2(n,1),Omega2(n,2))=tauu*D(Omega2(n,1), Omega2(n,2));
end


%%%%%%%%%%%%%%%%%%%Solve the SDR Problem%%%%%%%%%%%%%%%%%%%
cvx_begin quiet
      variable Q(L*N1,L*N1,K) hermitian;   %Variable for N x K beamforming matrix
      variable x(K,1) ;  %slack variables
      variable y(K,1) ;  %slack variables
      variable t(K,1);   %slack variables for precoder
      
minimize sum(t) %%%total transmit power
     subject to
     
     for k=1:K
       Q(:,:,k)==Wsolution_temp(:,k)*Wsolution_temp(:,k)';
     end
     
     for k=1:K        %%%%%%%%%Bernstein Type Inequality Constraint (1)
         temp=0; 
         CC_temp=(CC(:,k));  %N x 1 intantaneous CSI vector 
         DD_temp=diag(DD(:,k));  % N x N statistical CSI indicator matrix
         
         for j=1:K   %caculate the variable B_k   
           temp=temp+Q(:,:,j);  
         end
         B=(1/r+1)*Q(:,:,k)-temp;
         
         trace(DD_temp*B*DD_temp)+real(CC_temp'*Hermitian(B)*CC_temp)-sqrt(-2*log(epsilon))*x(k)+y(k)*log(epsilon)-delta^2>=0;
         %%%%%%%%%%%%%%Bernstein Type Inequality Constraint (2) %%%%%%%%
             norm(DD_temp*B*([DD_temp,sqrt(2)*CC_temp]), 'fro')<=x(k);
  
         y(k)*eye(L*N1)+DD_temp*B*DD_temp==semidefinite(L*N1);   
         y(k)>=0;
     end
     
     for l=1:length(A_set)    %%%Active RRHs: Transmit Power Constraints
          temp2=0;
         for k=1:K
         temp2=temp2+trace(Q(N1*(A_set(l)-1)+1:N1*A_set(l),N1*(A_set(l)-1)+1:N1*A_set(l),k));
         end
          temp2<=P(A_set(l));
     end
     
     for k=1:K
         Q(:,:,k)==semidefinite(L*N1);
         trace(Q(:,:,k))<=t(k);
         t(k)>=0;
     end
 cvx_end
     
     %Analyze result and prepare the output variables.
     if  strfind(cvx_status,'Solved') 
         feasible=true;
         Wsolution=Wsolution_temp;
     %%%%%%%%%%%%%%%%%%%%%%%%%
     else
         feasible=false;
         Wsolution=[];
     end