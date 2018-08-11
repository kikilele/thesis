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
%     去掉每个点自己与自己的连线，故为sizeIS(1)-n,
%     其中m点与n点之间的连线会算两遍，故(sizeIS(1)-n)/2,
%     那么(sizeIS(1)-n)/2在除以n就表示平均度了。
    avg=(sizeIS(1)-n)/2/n;
    fprintf('r=%.4f,平均度=%.4f\n',r,avg);
    degree=[degree,avg]; 
    distance=[distance,r];
    result=[degree;distance];
    
    
end