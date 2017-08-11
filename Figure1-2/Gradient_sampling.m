function [Gsolution] = Gradient_sampling(H_samples, W, M, K, L, N1, delta, r)
%INPUT:
%H_samples       = Kt*Nt x Kr matrix with row index for transmit and column
%                            index receiver antennas, M# samples
%W                    =# Beamforming matrix
%M                    =# Samples 
%K                     =# Users
%L                     =# RRHs
%N1                   =# Antennas
% t                     = Parameter for conservative approximation
%r                     = QoS requirements


%OUTPUT
%Gsolution       =KN x 1 vector with entries as the gradients

Gsolution=zeros(L*N1*K, 1); temp=zeros(L*N1*K,1);
for m=1:M
    S1=zeros(K,1); S2=0; S=0; temp2=0; temp_G=zeros(L*N1*K, 1); H=H_samples(:,:,m);
    for kk=1:K  %%%Calculate S_k, k=1,...K
         temp1=0; 
    for j=1:K
        if j~=kk
    temp1=temp1+(W(:,j)'*Hermitian(H(:,kk)*H(:,kk)')*W(:,j))+(1/r)*(W(:,j)'*Hermitian(H(:,j)*H(:,j)')*W(:,j));
        end
    end
    S1(kk)=temp1+delta^2;
    temp2=temp2+(1/r)*(W(:,kk)'*Hermitian(H(:,kk)*H(:,kk)')*W(:,kk));
    end
    %%%%Calculate S_(K+1)
    S2=temp2;
    S=[S1;S2]; % obtain the values S_k, k=1,...,K+1
    %%%%%%%%%% Calculate the gradient %%%%%%%%%%
    [Value, Index]=max(S);
    if Index<(K+1)
        
        for j=1:K
            if j~=Index
            temp_G(L*N1*(j-1)+1: L*N1*j, 1)=H(:,Index)*H(:,Index)'*W(:,j)+(1/r)*H(:,j)*H(:,j)'*W(:,j);
            end
        end
        
      else
        
        for j=1:K
        temp_G(L*N1*(j-1)+1: L*N1*j, 1)=(1/r)*H(:,j)*H(:,j)'*W(:,j);
        end
        
    end
   temp=temp+temp_G;   
end
 Gsolution=(1/M)*temp;
    