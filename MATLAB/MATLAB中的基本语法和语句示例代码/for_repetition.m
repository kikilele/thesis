clear
sum=zeros(6,1); %����0����6��1��
for i=eye(6,6)     %ÿ��ѭ��n����ȡeye(6,6)��һ��
    sum=sum+i;
end
fprintf('sum is %d\n',sum);
