function P = dgraf(A)
n = size(A,1); %����A������
P = A;

for i=2:n
    P = P + A^i;
end

P(P~=0) = 1;  %��������0��ֵ�޸�Ϊ1
P;