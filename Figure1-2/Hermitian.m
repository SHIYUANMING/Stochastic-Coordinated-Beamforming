function [H]=Hermitian(A)
D=size(A,1); H=A;
for i=1:D
    for j=i:D
        H(j,i)=A(i,j)';
    end
end
        