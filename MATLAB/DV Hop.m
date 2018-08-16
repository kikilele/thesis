%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  DV-Hop算法  ~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% BorderLength-----正方形区域的边长，单位：m
% NodeAmount-------网络节点的个数
% BeaconAmount---信标节点数
% Sxy--------------用于存储节点的序号，横坐标，纵坐标的矩阵
%Beacon----------信标节点坐标矩阵;BeaconAmount*BeaconAmount
%UN-------------未知节点坐标矩阵;2*UNAmount
% Distance------未知节点到信标节点距离矩阵;2*BeaconAmount
%h---------------节点间初始跳数矩阵
%X---------------节点估计坐标初始矩阵,X=[x,y]'
% R------------------节点的通信距离，一般为10-100m

clear,close all;
BorderLength=100;
NodeAmount=100;
BeaconAmount=8;
UNAmount=NodeAmount-BeaconAmount;
R=50;
% D=zeros(NodeAmount,NodeAmount);%未知节电到信标节点距离初始矩阵；BeaconAmount行NodeAmount列
h=zeros(NodeAmount,NodeAmount);%初始跳数为0；BeaconAmount行NodeAmount列
X=zeros(2,UNAmount);%节点估计坐标初始矩阵

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~在正方形区域内产生均匀分布的随机拓扑~~~~~~~~~~~~~~~~~~~~
C=BorderLength.*rand(2,NodeAmount);
%带逻辑号的节点坐标
Sxy=[[1:NodeAmount];C];
Beacon=[Sxy(2,1:BeaconAmount);Sxy(3,1:BeaconAmount)];%信标节点坐标
UN=[Sxy(2,(BeaconAmount+1):NodeAmount);Sxy(3,(BeaconAmount+1):NodeAmount)];%未知节点坐标
%画出节点分布图
plot(Sxy(2,1:BeaconAmount),Sxy(3,1:BeaconAmount),'r*',Sxy(2,(BeaconAmount+1):NodeAmount),Sxy(3,(BeaconAmount+1):NodeAmount),'k.')
xlim([0,BorderLength]);
ylim([0,BorderLength]);
title('* 红色信标节点 . 黑色未知节点')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~初始化节点间距离、跳数矩阵~~~~~~~~~~~~~~~~~~~~~~
for i=1:NodeAmount
    for j=1:NodeAmount
        Dall(i,j)=((Sxy(2,i)-Sxy(2,j))^2+(Sxy(3,i)-Sxy(3,j))^2)^0.5;%所有节点间相互距离
        if (Dall(i,j)<=R)&(Dall(i,j)>0)
            h(i,j)=1;%初始跳数矩阵
        elseif i==j
            h(i,j)=0;
        else h(i,j)=inf;
        end
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~最短路经算法计算节点间跳数~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k=1:NodeAmount
    for i=1:NodeAmount
        for j=1:NodeAmount
            if h(i,k)+h(k,j)<h(i,j)%min(h(i,j),h(i,k)+h(k,j))
                h(i,j)=h(i,k)+h(k,j);
            end
        end
    end
end
h
%~~~~~~~~~~~~~~~~~~~~~~~~~求每个信标节点的校正值~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h1=h(1:BeaconAmount,1:BeaconAmount); 
D1=Dall(1:BeaconAmount,1:BeaconAmount);
for i=1:BeaconAmount
    dhop(i,1)=sum(D1(i,:))/sum(h1(i,:));%每个信标节点的平均每跳距离
end
D2=Dall(1:BeaconAmount,(BeaconAmount+1):NodeAmount);%BeaconAmount行UNAmount列
for i=1:BeaconAmount
    for j=1:UNAmount
        if min(D2(:,j))==D2(i,j)
            Dhop(1,j)=D2(i,j);%未知节点从最近的信标获得校正值
        end
    end
end
Dhop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~用跳数估计距离~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hop1=h(1:BeaconAmount,(BeaconAmount+1):NodeAmount)%未知节点到信标跳数，BeaconAmount行UNAmount列
for i=1:UNAmount
    hop=Dhop(1,i);%hop为从最近信标获得的校正值
    Distance(:,i)=hop*hop1(:,i);%%Beacon行UN列；
end
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~最小二乘法求未知点坐标~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d=Distance;
for i=1:2
    for j=1:(BeaconAmount-1)
      a(i,j)=Beacon(i,j)-Beacon(i,BeaconAmount);
    end
end
A=-2*(a');
% d=d1';
 for m=1:UNAmount 
     for i=1:(BeaconAmount-1)
         B(i,1)=d(i,m)^2-d(BeaconAmount,m)^2-Beacon(1,i)^2+Beacon(1,BeaconAmount)^2-Beacon(2,i)^2+Beacon(2,BeaconAmount)^2;
     end
           X1=inv(A'*A)*A'*B;
           X(1,m)=X1(1,1);
           X(2,m)=X1(2,1);
 end
 UN
 X
 for i=1:UNAmount
     error(1,i)=(((X(1,i)-UN(1,i))^2+(X(2,i)-UN(2,i))^2)^0.5);
 end
 figure;plot(error,'-o')
 title('每个未知节点的误差')
 error=sum(error)/UNAmount
 Accuracy=error/R