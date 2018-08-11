clear
sum=zeros(6,1); %产生0矩阵，6行1列
for i=eye(6,6)     %每次循环n依次取eye(6,6)的一列
    sum=sum+i;
end
fprintf('sum is %d\n',sum);
