function P = dgraf(A)
n = size(A,1); %返回A的行数
P = A;

for i=2:n
    P = P + A^i;
end

P(P~=0) = 1;  %将不等于0的值修改为1
P;