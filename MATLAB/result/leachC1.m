
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ƽ������  
xm=100;
ym=100;

%��վ���ݺ�����
sink.x=0.5*xm;
sink.y=0.5*ym;

%�ڵ�����
n=100;

%��ͷ����
p=0.05;

%���е�������
rmax=2999;

%����ģ�ͣ���λ����      
%ÿһ���ڵ�ĳ�ʼ���ڵ�
Eo=0.2;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%�Ŵ�������������
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%�ۼ�������Ҫ���ĵ�����
EDA=5*0.000000001;
%initial
A=zeros(1,rmax);
%����һ������,���ڴ��ÿһ�����нڵ�ʣ������֮��.
D=zeros(1,rmax);
%�����ڵ����

%���ڴ��ÿ���л�վ�յ������ݰ�
T=zeros(1,rmax);

%�洢����BS�������Զ�������ڵ��������һά����
MIN = zeros(1,rmax+1);
MAX = zeros(1,rmax+1);



%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

%������ĵĽ���ֵ�������������friss free space model��С���������two-ray ground model 
do=sqrt(Efs/Emp);

%��������ڵ�
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %��ʱ���нڵ㶼����ͨ�ڵ�
    S(i).type='N';
 
   
    %���ڵ㸳��ֵ
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
 
%��������
figure(1);

%��ͷ������
countCHs=0;
%ÿһ���ֵĴ�ͷ��
rcountCHs=0;
cluster=1;

countCHs;
rcountCHs=rcountCHs+countCHs;
flag_first_dead=0;

for r=1:1:rmax
    r
%������
E=0;
  %ȫ���ڵ㶼������ͷ�����һ���ֻ�
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end

hold off;

%��������
dead=0;
%
dead_a=0;
%
dead_n=0;
%���ڵ���
alive=n;

%���͵���ͷ�ͻ�վ�����ݰ�
%��λ����
packets_TO_BS=0;
packets_TO_CH=0;
%ÿ�ַ��͵���ͷ�ͻ�վ�����ݰ�
%��λ����
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
    %����Ƿ��нڵ�������
    if (S(i).E<=0) %�ڵ��Ѿ�����
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

%��һ���ڵ㿪ʼ����
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end

countCHs=0;
cluster=1;
% leachc,������������С
for i=1:n;
    if(S(i).E>0)
    E=E+S(i).E;
    end
end


for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)

 %��ͷѡ��
 if(temp_rand<= ((S(i).E*alive*p)/E))  %С������ֵ����ʱ��Ϊ��ͷ
            countCHs=countCHs+1;
%             packets_TO_BS=packets_TO_BS+1;
%             PACKETS_TO_BS(r+1)=packets_TO_BS;
            
            S(i).type='C';%����Ѿ�������ͷ��
            S(i).G=round(1/p)-1; 
            %����Ϊ��ͷ�ڵ�������¼�������Ա���������ͨ�ڵ�ѡ����õĴ�ͷ
            C(cluster).xd=S(i).xd;%������
            C(cluster).yd=S(i).yd;%������
            plot(S(i).xd,S(i).yd,'r*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%�ýڵ㵽BS�ľ���
            C(cluster).distance=distance;
            %�Ѵ�ͷID��¼����
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            %��ͷ����1
            cluster=cluster+1; 
            
            %�����������
            distance;
            if (distance>do) 
                %����Friss free space model
                S(i).E=S(i).E- ( (ETX+EDA)*(2000) + Emp*2000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                %����Two-Ray Ground model
                S(i).E=S(i).E- ( (ETX+EDA)*(2000)  + Efs*2000*( distance * distance )); 
            end
        end     
    
    end
  end 
end

STATISTICS(r+1).CLUSTERHEADS=cluster-1;% cluster-1��ֵΪ�����еĴ�ͷ����
CLUSTERHS(r+1)=cluster-1;

%ѡ������ĸ���ͷ
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 ) %�жϸýڵ㲻�Ǵ�ͷ����û������
     if(cluster-1>=1)  %��ͷ���ڵ���1��
     %  min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );%�ڵ㵽��վ֮��ľ���
       min_dis=10000;%����Ϊ���
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       
       %��ͨ�ڵ�����
            min_dis;  %����һ������֪���ýڵ����һ���ǻ�վ���Ǵ�ͷ������Ǵ�ͷ�Ļ� ��Ҫ�ô�ͷת����ȥ
         if(min_dis>0)
             if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(2000) + Emp*2000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(2000) + Efs*2000*( min_dis * min_dis)); 
            end                                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %��ͨ�ڵ㴫����ͷ�ڵ㣬��ͷ�ڵ��ٷ���BS
            dd=C(min_dis_cluster).distance;
            if(dd>do)
                %S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*2000 + Emp*2000*( dd*dd*dd*dd ) );
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*2000 + Emp*2000*( dd*dd*dd*dd ) );
            end
            if(dd<=do)
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E - ( (ERX + EDA)*2000 + Efs*2000*( dd * dd));
            end
         PACKETS_TO_CH(r+1)=n-dead-cluster+1;%��ͷ�յ������ݰ�
         PACKETS_TO_BS(r+1)=n-dead-cluster+1;%��վ�ܵ������ݰ�
        end

       packets_TO_BS=packets_TO_BS+1;    
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
   
   end
   %��ͷ��Ϊ0ʱ�� Ӧ�� 
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
%������������нڵ������֮��
for i=1:1:n
    if(S(i).E>=0)
    A(1,r)=A(1,r)+S(i).E;
    end
end
%����ÿ���е�����������ÿ���ڵ�һ�η���1�����ݰ�(2000����)
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