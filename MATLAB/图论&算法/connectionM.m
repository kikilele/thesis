%这个函数可以产生任意点数的随机拓扑图以及对应的邻接矩阵
%参数为num，点数

function [L]=connectionM(num)
%产生数组A用来存放表示两点之间权值的矩阵A，也就是邻接矩阵，那么两点之间权值不为零元素的个数即为该点的度数
DEF=num;
L=rand(DEF);%随机产生一个DEF阶矩阵，然后以此为根据修改为邻接矩阵
for i=1:DEF
    L(i,i)=0;%对角线上元素改为0，无向图的邻接矩阵的性质
end
L=10*L;
L=floor(L);%下取整
L=mod(L,2);%邻接矩阵中的元素取值为0，1，2
for i=1:DEF
    for j=1:i
        L(j,i)=L(i,j);
    end
end
x=100*rand(1,DEF);
y=100*rand(1,DEF);
plot(x,y,'r+');
for i=1:DEF
    a=find(L(i,:)>0);%将L矩阵每行大于0的数的在该行的地址找出来放在a中
    for j=1:length(a)        
        hold on;
        line([x(i),x(a(j))],[y(i),y(a(j))]);
    end
end
title('随机拓扑图');
e=num2str(DEF); 
legend(e);
for m=1:DEF
    f=num2str(m); 
    hold on;
    text((x(m)+x(m))/2,(y(m)+y(m))/2,f,'Fontsize',18); 
end
hold off