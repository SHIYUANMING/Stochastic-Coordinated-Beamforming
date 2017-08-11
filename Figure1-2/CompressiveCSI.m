function [Omega1, Omega2] = CompressiveCSI(D, M)
%INPUT
%D:  Large-Scale Coefficients
%M: # Channel Coefficients That Can not be Obtained
K=size(D,2);
N=size(D,1);
M1=M;
M2=N-M1;
for k=1:K
    temp=D(:,k);
    sort_value=sort(temp);
    
    min_value=sort_value(1:M1);
    for m=1:M1
    [row,column]=find(D==min_value(m));   
    Omega1((k-1)*M1+m,1)=row;
    Omega1((k-1)*M1+m,2)=column;
    end 
    
    min_value=sort_value(M1+1:N);
    for m=1:M2
    [row,column]=find(D==min_value(m));   
    Omega2((k-1)*M2+m,1)=row;
    Omega2((k-1)*M2+m,2)=column;
    end 
    
 
end

