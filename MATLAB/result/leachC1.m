
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%平面坐标  
xm=100;
ym=100;

%基站的纵横坐标
sink.x=0.5*xm;
sink.y=0.5*ym;

%节点总数
n=100;

%簇头比例
p=0.05;

%运行的总轮数
rmax=2999;

%能力模型：单位焦耳      
%每一个节点的初始化节点
Eo=0.2;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%放大器的两个参数
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%聚集数据所要消耗的能量
EDA=5*0.000000001;
%initial
A=zeros(1,rmax);
%定义一个矩阵,用于存放每一轮所有节点剩余能量之和.
D=zeros(1,rmax);
%死亡节点个数

%用于存放每轮中基站收到的数据包
T=zeros(1,rmax);

%存储距离BS最近和最远的两个节点的能量的一维数组
MIN = zeros(1,rmax+1);
MAX = zeros(1,rmax+1);



%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%能量损耗的界限值：大于它则符合friss free space model，小于它则符合two-ray ground model 
do=sqrt(Efs/Emp);

%产生随机节点
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %此时所有节点都是普通节点
    S(i).type='N';
 
   
    %给节点赋初值
        S(i).E=Eo;
        S(i).ENERGY=0;
        plot(S(i).xd,S(i).yd,'o');
        hold on;
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'x');
    
 mind=sqrt( (S(1).xd-(S(n+1).xd) )^2 + (S(1).yd-(S(n+1).yd) )^2 );
 maxd=sqrt( (S(1).xd-(S(n+1).xd) )^2 + (S(1).yd-(S(n+1).yd) )^2 );
 min_id=1;
 max_id=1;
 for i=2:1:n
    tempd=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
    if(tempd>maxd)
    maxd=tempd;
    max_id=i;
    end
    if(tempd<mind)
    mind=tempd;
    min_id=i;
    end
 end
 
%产生窗口
figure(1);

%簇头计数器
countCHs=0;
%每一个轮的簇头数
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=1:1:rmax
    r
%总能量
E=0;
  %全部节点都作过簇头后的又一个轮回
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

hold off;

%死亡节数
dead=0;
%
dead_a=0;
%
dead_n=0;
%存活节点数
alive=n;

%发送到簇头和基站的数据包
%单位比特
packets_TO_BS=0;
packets_TO_CH=0;
%每轮发送到簇头和基站的数据包
%单位比特
PACKETS_TO_CH(r+1)=0;
PACKETS_TO_BS(r+1)=0;

figure(2);

if(S(min_id).E>0)
   MIN(r+1)=S(min_id).E; 
end
if(S(max_id).E>0)
   MAX(r+1)=S(max_id).E;
end

for i=1:1:n
    %检查是否有节点死亡了
    if (S(i).E<=0) %节点已经死亡
        S(i).E=0;
        plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        alive=alive-1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    end
    if S(i).E>0
        S(i).type='N';
        if (S(i).ENERGY==0)  
        plot(S(i).xd,S(i).yd,'o');
        end
        if (S(i).ENERGY==1)  
        plot(S(i).xd,S(i).yd,'+');
        end
       hold on;
    end
end
plot(S(n+1).xd,S(n+1).yd,'rx');


STATISTICS(r).DEAD=dead;
DEAD(r)=dead;
D(r)=dead;
DEAD_N(r)=dead_n;
DEAD_A(r)=dead_a;
ALIVE(r)=alive;

%第一个节点开始死亡
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
% leachc,计算总能量大小
for i=1:n;
    if(S(i).E>0)
    E=E+S(i).E;
    end
end


for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %簇头选举
 if(temp_rand<= ((S(i).E*alive*p)/E))  %小于门限值，此时作为簇头
            countCHs=countCHs+1;
%             packets_TO_BS=packets_TO_BS+1;
%             PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';%标记已经作过簇头了
            S(i).G=round(1/p)-1; 
            %把作为簇头节点的坐标记录下来。以便于其他普通节点选择最好的簇头
            C(cluster).xd=S(i).xd;%横坐标
            C(cluster).yd=S(i).yd;%纵坐标
            plot(S(i).xd,S(i).yd,'r*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%该节点到BS的距离
            C(cluster).distance=distance;
            %把簇头ID记录下来
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            %簇头数加1
            cluster=cluster+1; 
            
            %计算能量损耗
            distance;
            if (distance>do) 
                %符合Friss free space model
                S(i).E=S(i).E- ( (ETX+EDA)*(2000) + Emp*2000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                %符合Two-Ray Ground model
                S(i).E=S(i).E- ( (ETX+EDA)*(2000)  + Efs*2000*( distance * distance )); 
            end
        end     
    
    end
  end 
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;% cluster-1的值为该轮中的簇头数量
CLUSTERHS(r+1)=cluster-1;

%选择加入哪个簇头
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 ) %判断该节点不是簇头并且没有死亡
     if(cluster-1>=1)  %簇头大于等于1个
     %  min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%节点到基站之间的距离
       min_dis=10000;%设置为最大
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %普通节点的损耗
            min_dis;  %到这一步并不知道该节点的下一跳是基站还是簇头，如果是簇头的话 还要让簇头转发出去
         if(min_dis>0)
             if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(2000) + Emp*2000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(2000) + Efs*2000*( min_dis * min_dis)); 
            end                                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %普通节点传至簇头节点，簇头节点再发至BS
            dd=C(min_dis_cluster).distance;
            if(dd>do)
                %S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*2000 + Emp*2000*( dd*dd*dd*dd ) );
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*2000 + Emp*2000*( dd*dd*dd*dd ) );
            end
            if(dd<=do)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*2000 + Efs*2000*( dd * dd));
            end
         PACKETS_TO_CH(r+1)=n-dead-cluster+1;%簇头收到的数据包
         PACKETS_TO_BS(r+1)=n-dead-cluster+1;%基站受到的数据包
        end

       packets_TO_BS=packets_TO_BS+1;    
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
   
   end
   %簇头数为0时候 应该 
   if(cluster-1==0)
      ddd0 = sqrt((S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2);
            if (ddd0>do)
                S(i).E=S(i).E- ( ETX*(2000) + Emp*2000*( ddd0 * ddd0 * ddd0 * ddd0)); 
            end
            if (ddd0<=do)
                S(i).E=S(i).E- ( ETX*(2000) + Efs*2000*( ddd0 * ddd0)); 
            end
            packets_TO_BS=packets_TO_BS+1;
   end
 end
end
hold on;
%计算该轮中所有节点的能量之和
for i=1:1:n
    if(S(i).E>=0)
    A(1,r)=A(1,r)+S(i).E;
    end
end
%计算每轮中的总吞吐量。每个节点一次发送1个数据包(2000比特)
T(1,r)=T(1,r)+packets_TO_BS;

countCHs;
rcountCHs=rcountCHs+countCHs;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   STATISTICS    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                                     %
%  DEAD  : a rmax x 1 array of number of dead nodes/round 
%  DEAD_A : a rmax x 1 array of number of dead Advanced nodes/round
%  DEAD_N : a rmax x 1 array of number of dead Normal nodes/round
%  CLUSTERHS : a rmax x 1 array of number of Cluster Heads/round
%  PACKETS_TO_BS : a rmax x 1 array of number packets send to Base Station/round
%  PACKETS_TO_CH : a rmax x 1 array of number of packets send to ClusterHeads/round
%  first_dead: the round where the first node died                   
%                                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%