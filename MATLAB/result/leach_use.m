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
packetLength = 4000;%数据包长度
ctrPacketLength = 100;%控制包长度
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
    S(i).xd=rand(1,1)*xm;%坐标
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';%普通节点
  
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

for r=0:1:rmax %主循环,每次1轮
r;

%Operation for epoch
if(mod(r, round(1/p2) )==0)%round(1/p)取整，四舍五入
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
if (dead2 == n)%节点全部死亡退出循环
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
     if ( (S(i).G)<=0) %如果该节点在候选集合中
        %Election of Cluster Heads
        if( temp_rand2 <= (p2/(1-p2*mod(r,round(1/p2)))))
            countCHs2 = countCHs2+1;
           
            S(i).type = 'C';
            S(i).G = round(1/p2)-1;
            CL2(cluster2).xd = S(i).xd;
            CL2(cluster2).yd = S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
           
            distance2=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%到sink的距离
            CL2(cluster2).distance = distance2;
            CL2(cluster2).id = i;
            X1(cluster2)=S(i).xd;
            Y1(cluster2)=S(i).yd;
            cluster2=cluster2+1;
            %%广播自成为簇头
            distanceBroad2 = sqrt(xm*xm+ym*ym);
            if (distanceBroad2 > do)
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Emp* ctrPacketLength*( distanceBroad2*distanceBroad2*distanceBroad2*distanceBroad2 ));%广播自成为簇头
            else
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Efs * ctrPacketLength*( distanceBroad2*distanceBroad2));
            end
            %Calculation of Energy dissipated 簇头自己发送数据包能量消耗
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

STATISTICS2(r+1).CLUSTERHEADS = cluster2-1;%统计第r轮簇头数目,r是从0开始的,所以加1;cluster最后要-1,是应为上面的循环多加了1
CLUSTERHS2(r+1)= cluster2-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 ) %普通节点
    % min_dis = sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%默认距离是到sink的距离
     min_dis2 = INFINITY;
     if(cluster2 -1 >= 1)%如果有簇头存在
         min_dis_cluster2 = 1;
         %加入最近的簇头
         for c = 1:1:cluster2 - 1 %簇头数量一共是cluster - 1
            %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
            temp2 = sqrt( (S(i).xd - CL2(c).xd)^2 + (S(i).yd - CL2(c).yd)^2 );
            if ( temp2 < min_dis2 )
                min_dis2 = temp2;
                min_dis_cluster2 = c;
            end
            %接收簇头发来的广播的消耗
            S(i).E = S(i).E - ETX * ctrPacketLength;
         end
      
         %Energy dissipated by associated Cluster Head普通节点发送数据包到簇头消耗,和加入消息
         min_dis2;
         if (min_dis2 > do)
             S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2)); %向簇头发送加入控制消息
             S(i).E = S(i).E - ( ETX*(packetLength) + Emp*packetLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2)); %向簇头数据包
         else
            S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis2 * min_dis2)); %向簇头发送加入控制消息
            S(i).E = S(i).E - ( ETX*(packetLength) + Efs*packetLength*( min_dis2 * min_dis2)); %向簇头数据包
         end
         S(i).E = S(i).E - ETX*(ctrPacketLength); %接收簇头确认加入控制消息
            
         %Energy dissipated %簇头接收簇成员数据包消耗能量,接收加入消息和和确认加入消息
         if(min_dis2 > 0)
            S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( (ERX + EDA)*packetLength ); %接受簇成员发来的数据包
            S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ERX *ctrPacketLength ; %接收加入消息
            if (min_dis2 > do)%簇头向簇成员发送确认加入的消息
                S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2));
            else
                S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis2 * min_dis2));
            end
           PACKETS_TO_CH2(r+1) = n - dead2 - cluster2 + 1; %所有的非死亡的普通节点都发送数据包
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
title('图一：存活节点数统计图');
xlabel('工作轮数（/轮）');
ylabel('存活节点数（/个）');
hold off
%r=0:1:rmax;
%plot(r,EJ(r+1));
%box off;
%title('图三：网络剩余能量统计图');
%xlabel('工作轮数（/轮）');
%ylabel('网络剩余能量（/焦耳）');

