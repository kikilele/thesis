clear;
clc;
nodenumber=1000;  %100��nodes
length=nodenumber;
width=nodenumber;
sink_x=nodenumber;  %sink�ڵ��λ��(100,100)
sink_y=nodenumber;
source_x=0; %Դ�ڵ��λ��(0,0)
source_y=0;
sinkposition=[length;width];
energyinitialvalue=1;
circle=10;
reliabletransmissiondistance=30;
packetnumber=21;  %���ݰ�����
packetlength=4000; % ���ݰ���С Data packet size
% delay��ȷ��
Ttr=0.4e-003;  %����ÿ�ֽ�����Ҫ��ʱ��
Tre=0.4e-003;  %����ÿ�ֽ�����Ҫ��ʱ��
Delay=50*(Ttr+Tre);  %�������ݰ�ʱ���ӳ�
%Radio amplifer energy;
Eelec=50*0.000000001;
Efs=10*0.000000000001;  
Eamp=0.0013*0.000000000001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ڵ�λ����Ϣ
%������Դ�ڵ��sink�ڵ����ڵĽڵ�λ�þ���
nodeposition(1,1)=source_x; 
nodeposition(2,1)=source_y;
nodeposition(1,nodenumber)=sink_x; 
nodeposition(2,nodenumber)=sink_y;
randlocation=nodenumber*rand(2,nodenumber);  %�������ڵ�
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