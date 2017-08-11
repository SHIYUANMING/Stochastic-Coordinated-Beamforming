function [H_samples] = samples(H, D, Omega1, Omega2, N1, T, tauu)
%INPUT:
%Omega:  |Omega| x 2 matrix Positions of CSI that can only obtain statistical CSI
%N1:       #antennas at each RRH
%T:          # Samples
M1=size(Omega1,1); 
M2=size(Omega2,1);

for tt=1:T   
    
for n=1:M1
H_temp(Omega1(n,1),Omega1(n,2))=D(Omega1(n,1), Omega1(n,2))*(normrnd(0,1/sqrt(2),N1,1)+i*normrnd(0,1/sqrt(2),N1,1));
end
H_samples(:,:,tt)=H_temp;

for n=1:M2
H_temp(Omega2(n,1),Omega2(n,2))=sqrt(1-tauu^2)*H(Omega2(n,1), Omega2(n,2))+D(Omega2(n,1), Omega2(n,2))*(tauu*(normrnd(0,1/sqrt(2),N1,1)+i*normrnd(0,1/sqrt(2),N1,1)));
end
H_samples(:,:,tt)=H_temp;
end