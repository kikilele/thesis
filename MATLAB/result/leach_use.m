clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm = 30;
ym = 30;

%x and y Coordinates of the Sink
sink.x =0.5 * xm;
sink.y =0.5 * ym;

%Number of Nodes in the field
n = 30;

%Optimal Election Probability of a node
%to become cluster head
p2=0.05;
packetLength = 4000;%���ݰ�����
ctrPacketLength = 100;%���ư�����
%Energy Model (all values in Joules)
%Initial Energy
Eo = 0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;

%Values for Hetereogeneity
%Percentage of nodes than are advanced
m=1;
%\alpha
a=1;
INFINITY = 999999999999999;
%maximum number of rounds
rmax=2;

%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%Computation of do
do=sqrt(Efs/Emp);

%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;%����
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';%��ͨ�ڵ�
  
    %temp_rnd0 = i;
    %Random Election of Normal Nodes
    %if (temp_rnd0>=m*n+1)
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
   % end
    %Random Election of Advanced Nodes
   % if (temp_rnd0<m*n+1)
   %     S(i).E=Eo*(1+a)
   %     S(i).ENERGY=1;
   %     %plot(S(i).xd,S(i).yd,'+');
    %    %hold on;
   % end
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
   
       
%First Iteration
figure(1);

%counter for CHs
countCHs2=0;
%counter for CHs per round
rcountCHs2=0;
cluster2=1;

countCHs2;
rcountCHs2=rcountCHs2+countCHs2;
flag_first_dead2=0;

for r=0:1:rmax %��ѭ��,ÿ��1��
r;

%Operation for epoch
if(mod(r, round(1/p2) )==0)%round(1/p)ȡ������������
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
end

hold off;

%Number of dead nodes
dead2=0;
%Number of dead Advanced Nodes
dead_a2=0;
%Number of dead Normal Nodes
dead_n2=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS2=0;
packets_TO_CH2=0;
%counter for bit transmitted to Bases Station and to Cluster Heads
%per round
PACKETS_TO_CH2(r+1)=0;
PACKETS_TO_BS2(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
       % plot(S(i).xd,S(i).yd,'red .');
        dead2=dead2+1;
     end
    if S(i).E>0
        S(i).type='N';
    
    end
end
%plot(S(n+1).xd,S(n+1).yd,'x');
if (dead2 == n)%�ڵ�ȫ�������˳�ѭ��
   break;
end

STATISTICS2(r+1).DEAD=dead2;
DEAD2(r+1)=dead2;
DEAD_N2(r+1)=dead_n2;
DEAD_A2(r+1)=dead_a2;

%When the first node dies
if (dead2==1)
    if(flag_first_dead2==0)
        first_dead2=r;
        flag_first_dead2=1;
    end
end


countCHs2=0;
cluster2=1;
for i=1:1:n
   if(S(i).E>0)
     temp_rand2=rand;    
     if ( (S(i).G)<=0) %����ýڵ��ں�ѡ������
        %Election of Cluster Heads
        if( temp_rand2 <= (p2/(1-p2*mod(r,round(1/p2)))))
            countCHs2 = countCHs2+1;
           
            S(i).type = 'C';
            S(i).G = round(1/p2)-1;
            CL2(cluster2).xd = S(i).xd;
            CL2(cluster2).yd = S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
           
            distance2=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%��sink�ľ���
            CL2(cluster2).distance = distance2;
            CL2(cluster2).id = i;
            X1(cluster2)=S(i).xd;
            Y1(cluster2)=S(i).yd;
            cluster2=cluster2+1;
            %%�㲥�Գ�Ϊ��ͷ
            distanceBroad2 = sqrt(xm*xm+ym*ym);
            if (distanceBroad2 > do)
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Emp* ctrPacketLength*( distanceBroad2*distanceBroad2*distanceBroad2*distanceBroad2 ));%�㲥�Գ�Ϊ��ͷ
            else
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Efs * ctrPacketLength*( distanceBroad2*distanceBroad2));
            end
            %Calculation of Energy dissipated ��ͷ�Լ��������ݰ���������
            distance2;
            if (distance2 > do)
                 S(i).E = S(i).E- ( (ETX+EDA)*(packetLength) + Emp * packetLength*( distance2*distance2*distance2*distance2 ));
            else
                 S(i).E = S(i).E- ( (ETX+EDA)*(packetLength) + Efs * packetLength*( distance2 * distance2 ));
            end
            packets_TO_BS2 = packets_TO_BS2+1;
            PACKETS_TO_BS2(r+1) = packets_TO_BS2;
        end    
   
     end
   end
end

STATISTICS2(r+1).CLUSTERHEADS = cluster2-1;%ͳ�Ƶ�r�ִ�ͷ��Ŀ,r�Ǵ�0��ʼ��,���Լ�1;cluster���Ҫ-1,��ӦΪ�����ѭ�������1
CLUSTERHS2(r+1)= cluster2-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 ) %��ͨ�ڵ�
    % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%Ĭ�Ͼ����ǵ�sink�ľ���
     min_dis2 = INFINITY;
     if(cluster2 -1 >= 1)%����д�ͷ����
         min_dis_cluster2 = 1;
         %��������Ĵ�ͷ
         for c = 1:1:cluster2 - 1 %��ͷ����һ����cluster - 1
            %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
            temp2 = sqrt( (S(i).xd - CL2(c).xd)^2 + (S(i).yd - CL2(c).yd)^2 );
            if ( temp2 < min_dis2 )
                min_dis2 = temp2;
                min_dis_cluster2 = c;
            end
            %���մ�ͷ�����Ĺ㲥������
            S(i).E = S(i).E - ETX * ctrPacketLength;
         end
      
         %Energy dissipated by associated Cluster Head��ͨ�ڵ㷢�����ݰ�����ͷ����,�ͼ�����Ϣ
         min_dis2;
         if (min_dis2 > do)
             S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2)); %���ͷ���ͼ��������Ϣ
             S(i).E = S(i).E - ( ETX*(packetLength) + Emp*packetLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2)); %���ͷ���ݰ�
         else
            S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis2 * min_dis2)); %���ͷ���ͼ��������Ϣ
            S(i).E = S(i).E - ( ETX*(packetLength) + Efs*packetLength*( min_dis2 * min_dis2)); %���ͷ���ݰ�
         end
         S(i).E = S(i).E - ETX*(ctrPacketLength); %���մ�ͷȷ�ϼ��������Ϣ
            
         %Energy dissipated %��ͷ���մس�Ա���ݰ���������,���ռ�����Ϣ�ͺ�ȷ�ϼ�����Ϣ
         if(min_dis2 > 0)
            S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( (ERX + EDA)*packetLength ); %���ܴس�Ա���������ݰ�
            S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ERX *ctrPacketLength ; %���ռ�����Ϣ
            if (min_dis2 > do)%��ͷ��س�Ա����ȷ�ϼ������Ϣ
                S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2));
            else
                S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis2 * min_dis2));
            end
           PACKETS_TO_CH2(r+1) = n - dead2 - cluster2 + 1; %���еķ���������ͨ�ڵ㶼�������ݰ�
         end
      
         S(i).min_dis = min_dis2;
         S(i).min_dis_cluster = min_dis_cluster2;
    
     end
end
end
%hold on;

countCHs2;
rcountCHs2 = rcountCHs2 + countCHs2;



%Code for Voronoi Cells
%Unfortynately if there is a small
%number of cells, Matlab's voronoi
%procedure has some problems

%[vx,vy]=voronoi(X,Y);
%plot(X,Y,'r*',vx,vy,'b-');
% hold on;
% voronoi(X,Y);
% axis([0 xm 0 ym]);

end
figure(2);
x=1:1:r;
y=1:1:r;
z=1:1:r;

for i=1:r;
    x(i)=i;
    y(i) = n - STATISTICS2(i).DEAD;
    z(i)=CLUSTERHS2(i);
end
%plot(x,y,'r',x,z,'b');
plot(x,y,'r--');

hold on;
title('ͼһ�����ڵ���ͳ��ͼ');
xlabel('����������/�֣�');
ylabel('���ڵ�����/����');
hold off
%r=0:1:rmax;
%plot(r,EJ(r+1));
%box off;
%title('ͼ��������ʣ������ͳ��ͼ');
%xlabel('����������/�֣�');
%ylabel('����ʣ��������/������');

