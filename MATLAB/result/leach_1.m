clear;
%1.初始参数设定模块
%.传感器节点区域界限(单位 M)
xm=30;
ym=30;
%(1)汇聚节坐标给定
sink.x=xm;
sink.y=ym;
%区域内传器节数
n=30;
%簇头优化比例（当选簇头的概率）
p=0.3;
%能量模型（单位 焦）
%初始化能量模型
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
%最大循环次数
rmax=1;
%算出参数 do
do=sqrt(Efs/Emp);
%2.无线传感器网络模型产生模块
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    S(i).E=Eo;
    %initially there are no cluster heads only nodes
    S(i).type='N';
    plot(S(i).xd,S(i).yd,'bo');
    hold on;
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'rx');
figure(1);
%3.网络运行模块
%簇头节点数
countCHs1=0;
cluster1=1;%此定义的目的仅仅是给定一个1开始的下标参数，真正的簇头数应该还减去1
flag_first_dead1=0;
flag_all_dead1=0;
%死亡节点数
dead1=0;
first_dead1=0;
all_dead1=0;
%活动节点数
allive1=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS1=0;
packets_TO_CH1=0;
%(1)循环模式设定
for r=0:1:rmax     %该 for 循环将下面的所有程序包括在内，直到最后一 end 才结束循环    
  %每过一个轮转周期使各节点的S（i）.G参数（该参数用于后面的簇选举，在该轮转周期内已当选过簇头的节点不能再当选）恢复为零
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end
%(2)死亡节点检查模块
dead=0;
for i=1:1:n
    %检查有无死亡节点
    if (S(i).E<=0)
        dead1=dead1+1; 
        plot(S(i).xd,S(i).yd,'g+');
        %(3)第一个死亡节点的产生时间(用轮次表示)
        %第一个节点死亡时间
        if (dead==1)
           if(flag_first_dead==0)
              first_dead=r;
              flag_first_dead=1;
           end
        end
        
        if(dead==n)
           if(flag_all_dead==0)
              all_dead=r;
              flag_all_dead=1;
           end
        end
    end
    if S(i).E>0
        S(i).type='N';
    end
end
STATISTICS.DEAD(r+1)=dead;
%STATISTICS.ALLIVE(r+1)=allive-dead;
%(4)簇头选举模块
countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)  
       %簇头的选举，当选的簇头会把各种相关信存入下面程序所给定的变量中
        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            C(cluster).id=i;
            X(cluster)=S(i).xd;%CH的x坐标
            Y(cluster)=S(i).yd;%CH的y坐标
            cluster=cluster+1;
        end
   end
   end
end
STATISTICS.COUNTCHS(r+1)=countCHs;
P=[X S(n+1).xd];
Q=[Y S(n+1).yd];
%根据do计算权重矩阵
for i=1:cluster
    for j=1:cluster
        D(i,j)=sqrt((P(i)-P(j))^2 + (Q(i)-Q(j))^2 );
        if D(i,j)>do
            D(i,j)=inf;   
        elseif i==j 
            D(i,j)=0;
        else D(i,j)=1;
        end
    end
end

for i=1:cluster
    clear a;
    a=find(D(i,:)==1);
    if ~isempty(a)
        for j=1:length(a)           
            hold on;             
            line([P(i) P(a(j))],[Q(i) Q(a(j))]);
        end  
    end    
end         


%多路径转发
[d,index1,index2]=short(D);

%(5)簇内成员选择簇头模块(即簇的形成模块)
%簇内成员对簇头的选择（即簇的形成）算法
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=Inf;
       min_dis_cluster=0;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       %簇内节点（发送4000bit数据）能量消耗
       
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %簇头（接受和融合这一簇内节点4000bit数据）的能量消耗
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
%            packets_TO_CH=packets_TO_CH+1;    
 
        S(i).min_dis=min_dis;
        S(i).min_dis_cluster=min_dis_cluster;
    else
        min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
     end
  end
end
STATISTICS.PACKETS_TO_CH(r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(r+1)=packets_TO_BS;
end
first_dead;
all_dead;
STATISTICS.DEAD(r+1)
STATISTICS.ALLIVE(r+1)
STATISTICS.PACKETS_TO_CH(r+1)
STATISTICS.PACKETS_TO_BS(r+1)
STATISTICS.COUNTCHS(r+1)
r=0:1;
figure(2);
grid on;
plot(r,STATISTICS.DEAD);
xlabel('x(time)');
ylabel('y(dead)');



