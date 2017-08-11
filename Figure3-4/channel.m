clc; clear all;
L=5; K=3; N1=1; %%% #RRH:L, #user: K, #antenna: N1 
CC1=1000; CC2=1;

for cc1=1:CC1
%%%%%%%%%%%%%%%%Network Realization%%%%%%%%%%%%%%%%%%%%%%%%
 U_position=600.*(rand(2,K)-0.5);  %% user positions
 B_position=600.*(rand(2,L)-0.5);  %%RRH positions

 %%%%%Generate Large-Scale Fading%%%%%%%%%%%%%
for k=1:K
    for l=1:L
                   d=(norm(B_position(:,l)-U_position(:,k))+10);
                   D_temp(l,k)=4.4*10^(5)/(d^(1.88)*10^(normrnd(0,6.3)/20));
    end
end

load('D_temp.mat')

for cc2=1:CC2
%%%%%%Generate Small-Scale Fading%%%%%%%%%%%%%
for k=1:K
    for l=1:L
            H_temp(N1*(l-1)+1:N1*l,k)=D_temp(l,k)*(normrnd(0,1/sqrt(2),N1,1)+i*normrnd(0,1/sqrt(2),N1,1));  %%%nosie normalized to 1
    end
end

H(:,:,cc1)=H_temp;
D(:,:,cc1)=D_temp;

end
end

save('H.mat','H');
save('D.mat','D');