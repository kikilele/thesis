clear;
clc;
nodenumber=1000;  %100个nodes
length=nodenumber;
width=nodenumber;
sink_x=nodenumber;  %sink节点的位置(100,100)
sink_y=nodenumber;
source_x=0; %源节点的位置(0,0)
source_y=0;
sinkposition=[length;width];
energyinitialvalue=1;
circle=10;
reliabletransmissiondistance=30;
packetnumber=21;  %数据包个数
packetlength=4000; % 数据包大小 Data packet size
% delay的确定
Ttr=0.4e-003;  %发送每字节所需要的时间
Tre=0.4e-003;  %接受每字节所需要的时间
Delay=50*(Ttr+Tre);  %发送数据包时间延迟
%Radio amplifer energy;
Eelec=50*0.000000001;
Efs=10*0.000000000001;  
Eamp=0.0013*0.000000000001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%节点位置信息
%做包含源节点和sink节点在内的节点位置矩阵
nodeposition(1,1)=source_x; 
nodeposition(2,1)=source_y;
nodeposition(1,nodenumber)=sink_x; 
nodeposition(2,nodenumber)=sink_y;
randlocation=nodenumber*rand(2,nodenumber);  %随机部署节点
nodeenergyinitial(1,1)=100;
nodeenergyinitial(1,nodenumber)=100;
for i=2:1:nodenumber-1
    nodeposition(1,i)=randlocation(1,i-1);
    nodeposition(2,i)=randlocation(2,i-1);
    plot(nodeposition(1,i),nodeposition(2,i),'r.','markersize',10);
    nodeenergyinitial(1,i)=energyinitialvalue;
    hold on;
end
plot(nodeposition(1,1),nodeposition(2,1),'b^-','markersize',16);
plot(sinkposition(1,1),sinkposition(2,1),'b^-','markersize',16);