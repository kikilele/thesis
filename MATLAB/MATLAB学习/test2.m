clear;

n=200;
p=rand(n,1)+1i*rand(n,1);

A=repmat(p,1,n)-repmat(p.',n,1);
D=abs(A);
save figuredate.mat p A D;
degree=[];
distance=[];

for r=0.2:0.001:0.3
    [IS,JS]=find(D<r);   
    
    sizeIS=size(IS);
%     ȥ��ÿ�����Լ����Լ������ߣ���ΪsizeIS(1)-n,
%     ����m����n��֮������߻������飬��(sizeIS(1)-n)/2,
%     ��ô(sizeIS(1)-n)/2�ڳ���n�ͱ�ʾƽ�����ˡ�
    avg=(sizeIS(1)-n)/2/n;
    fprintf('r=%.4f,ƽ����=%.4f\n',r,avg);
    degree=[degree,avg]; 
    distance=[distance,r];
    result=[degree;distance];
    
    
end