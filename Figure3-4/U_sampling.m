function [Usolution] = U_sampling(H_samples, W, M, K, Mu, delta, r)
%INPUT:
%H_samples       = Kt*Nt x Kr matrix with row index for transmit and column
%                            index receiver antennas, M# samples
%W                    =# Beamforming matrix
%M                    =# Samples 
%K                     =# Users
% Mu                     = Parameter for conservative approximation
% r                     = QoS requirements

%OUTPUT
%Usolution       =1 x 1 the average sampling of U

Usolution=0; temp=0;
for m=1:M
    S1=zeros(K,1); S=0; S2=0; temp2=0; H=H_samples(:,:,m);
    for kk=1:K  %%%Calculate S_k, k=1,...K
         temp1=0; 
    for j=1:K
        if j~=kk
    temp1=temp1+(W(:,j)'*Hermitian(H(:,kk)*H(:,kk)')*W(:,j))+(1/r)*(W(:,j)'*Hermitian(H(:,j)*H(:,j)')*W(:,j));
        end
    end
    S1(kk)=temp1+delta^2+Mu;
    temp2=temp2+(1/r)*(W(:,kk)'*Hermitian(H(:,kk)*H(:,kk)')*W(:,kk));
    end
    %%%%%%%%%%%Calculate S_(K+1)%%%%%%%%%%%%
    S2=temp2;
    S=[S1;S2]; %obtain the values S_k, k=1,...,K+1
    temp=temp+max(S);
end
Usolution=(1/M)*temp;