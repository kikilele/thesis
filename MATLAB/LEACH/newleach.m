clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Field Dimensions - x and y maximum (in meters)
xm = 100;
ym = 100;

%x and y Coordinates of the Sink
sink.x =0.5 * xm;
sink.y =0.5 * ym;

%Number of Nodes in the field
n = 100;

%Optimal Election Probability of a node
%to become cluster head
p=0.05;
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
rmax=4000;

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
countCHs=0;
%counter for CHs per round
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=0:1:rmax %��ѭ��,ÿ��1��
r;

%Operation for epoch
if(mod(r, round(1/p) )==0)%round(1/p)ȡ������������
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
end

hold off;

%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;

%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%counter for bit transmitted to Bases Station and to Cluster Heads
%per round
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(1);

for i=1:1:n
    %checking if there is a dead node
    if (S(i).E<=0)
       % plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
     end
    if S(i).E>0
        S(i).type='N';
    
    end
end
%plot(S(n+1).xd,S(n+1).yd,'x');
if (dead == n)%�ڵ�ȫ�������˳�ѭ��
   break;
end

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r;
        flag_first_dead=1;
    end
end


countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
     temp_rand=rand;    
     if ( (S(i).G)<=0) %����ýڵ��ں�ѡ������
        %Election of Cluster Heads
        if( temp_rand <= (p/(1-p*mod(r,round(1/p)))))
            countCHs = countCHs+1;
           
            S(i).type = 'C';
            S(i).G = round(1/p)-1;
            C(cluster).xd = S(i).xd;
            C(cluster).yd = S(i).yd;
            %plot(S(i).xd,S(i).yd,'k*');
           
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%��sink�ľ���
            C(cluster).distance = distance;
            C(cluster).id = i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
            %%�㲥�Գ�Ϊ��ͷ
            distanceBroad = sqrt(xm*xm+ym*ym);
            if (distanceBroad > do)
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Emp* ctrPacketLength*( distanceBroad*distanceBroad*distanceBroad*distanceBroad ));%�㲥�Գ�Ϊ��ͷ
            else
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Efs * ctrPacketLength*( distanceBroad*distanceBroad));
            end
            %Calculation of Energy dissipated ��ͷ�Լ��������ݰ���������
            distance;
            if (distance > do)
                 S(i).E = S(i).E- ( (ETX+EDA)*(packetLength) + Emp * packetLength*( distance*distance*distance*distance ));
            else
                 S(i).E = S(i).E- ( (ETX+EDA)*(packetLength) + Efs * packetLength*( distance * distance ));
            end
            packets_TO_BS = packets_TO_BS+1;
            PACKETS_TO_BS(r+1) = packets_TO_BS;
        end    
   
     end
   end
end

STATISTICS(r+1).CLUSTERHEADS = cluster-1;%ͳ�Ƶ�r�ִ�ͷ��Ŀ,r�Ǵ�0��ʼ��,���Լ�1;cluster���Ҫ-1,��ӦΪ�����ѭ�������1
CLUSTERHS(r+1)= cluster-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 ) %��ͨ�ڵ�
    % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%Ĭ�Ͼ����ǵ�sink�ľ���
     min_dis = INFINITY;
     if(cluster -1 >= 1)%����д�ͷ����
         min_dis_cluster = 1;
         %��������Ĵ�ͷ
         for c = 1:1:cluster - 1 %��ͷ����һ����cluster - 1
            %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
            temp = sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 );
            if ( temp < min_dis )
                min_dis = temp;
                min_dis_cluster = c;
            end
            %���մ�ͷ�����Ĺ㲥������
            S(i).E = S(i).E - ETX * ctrPacketLength;
         end
      
         %Energy dissipated by associated Cluster Head��ͨ�ڵ㷢�����ݰ�����ͷ����,�ͼ�����Ϣ
         min_dis;
         if (min_dis > do)
             S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis)); %���ͷ���ͼ��������Ϣ
             S(i).E = S(i).E - ( ETX*(packetLength) + Emp*packetLength*( min_dis * min_dis * min_dis * min_dis)); %���ͷ���ݰ�
         else
            S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis * min_dis)); %���ͷ���ͼ��������Ϣ
            S(i).E = S(i).E - ( ETX*(packetLength) + Efs*packetLength*( min_dis * min_dis)); %���ͷ���ݰ�
         end
         S(i).E = S(i).E - ETX*(ctrPacketLength); %���մ�ͷȷ�ϼ��������Ϣ
            
         %Energy dissipated %��ͷ���մس�Ա���ݰ���������,���ռ�����Ϣ�ͺ�ȷ�ϼ�����Ϣ
         if(min_dis > 0)
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*packetLength ); %���ܴس�Ա���������ݰ�
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ERX *ctrPacketLength ; %���ռ�����Ϣ
            if (min_dis > do)%��ͷ��س�Ա����ȷ�ϼ������Ϣ
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis * min_dis * min_dis * min_dis));
            else
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis * min_dis));
            end
           PACKETS_TO_CH(r+1) = n - dead - cluster + 1; %���еķ���������ͨ�ڵ㶼�������ݰ�
         end
      
         S(i).min_dis = min_dis;
         S(i).min_dis_cluster = min_dis_cluster;
    
     end
end
end
%hold on;

countCHs;
rcountCHs = rcountCHs + countCHs;



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

x=1:1:r;
y=1:1:r;
z=1:1:r;

for i=1:r;
    x(i)=i;
    y(i) = n - STATISTICS(i).DEAD;
    z(i)=CLUSTERHS(i);
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

