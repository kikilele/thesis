clear
%1.��ʼ�����趨ģ��
%.�������ڵ��������(��λ M)
xm=100;
ym=100;
%(1)��۽��������
sink.x=0.5*xm;
sink.y=0.5*ym;
%�����ڴ�������
n=100;
%��ͷ�Ż���������ѡ��ͷ�ĸ��ʣ�
p=0.1;
%����ģ�ͣ���λ ����
%��ʼ������ģ��
Eo=0.5;
%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
%�������ڵ㳬��һ�ڵ������İٷֱ�
a=1;
%���ѭ������
rmax=2;
%������� do
do=sqrt(Efs/Emp);
%2.���ߴ���������ģ�Ͳ���ģ��
%�������ߴ���������,�������ھ���Ͷ��100���ڵ�,������ͼ��
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    S(i).E=Eo*(1+rand*a);
    %initially there are no cluster heads only nodes
    S(i).type='N';
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%3.��������ģ��
%��ͷ�ڵ���
countCHs=0;
cluster=1;%�˶����Ŀ�Ľ����Ǹ���һ��1��ʼ���±�����������Ĵ�ͷ��Ӧ�û���ȥ1
flag_first_dead=0;
flag_teenth_dead=0;
flag_all_dead=0;
%�����ڵ���
dead=0;
first_dead=0;
teenth_dead=0;
all_dead=0;
%��ڵ���
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
%(1)ѭ��ģʽ�趨
for r=0:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r
  %ÿ��һ����ת����(������Ϊ10��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣��ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ��
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        S(i).cl=0;
    end
  end
%(2)�����ڵ���ģ��
dead=0;
for i=1:1:n
    %������������ڵ�
    if (S(i).E<=0)
        dead=dead+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead==1)
           if(flag_first_dead==0)
              first_dead=r;
              flag_first_dead=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead==0.1*n)
           if(flag_teenth_dead==0)
              teenth_dead=r;
              flag_teenth_dead=1;
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
STATISTICS.ALLIVE(r+1)=allive-dead;
%(4)��ͷѡ��ģ��
countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)  
       %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
            packets_TO_BS=packets_TO_BS+1;
            PACKETS_TO_BS(r+1)=packets_TO_BS;
             S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
           distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            cluster=cluster+1;
           %�����ͷ����4000bit���ݵ���վ���������ģ�����Ӧ�����нڵ������ͷÿһ�ַ���4000bit���ݣ�
           distance;
            if (distance>do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
    end
    % S(i).G=S(i).G-1;  
   
 end 
end
STATISTICS.COUNTCHS(r+1)=countCHs;
%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
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
       %���ڽڵ㣨����4000bit���ݣ���������
       
            min_dis;
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %��ͷ�����ܺ��ں���һ���ڽڵ�4000bit���ݣ�����������
            S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
            packets_TO_CH=packets_TO_CH+1;    
 
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
first_dead
teenth_dead
all_dead
STATISTICS.DEAD(r+1)
STATISTICS.ALLIVE(r+1)
STATISTICS.PACKETS_TO_CH(r+1)
STATISTICS.PACKETS_TO_BS(r+1)
STATISTICS.COUNTCHS(r+1)
r=0:5000;
subplot(2,2,1);
plot(r,STATISTICS.DEAD);
subplot(2,2,2);
plot(r,STATISTICS.ALLIVE);
subplot(2,2,3);
plot(r,STATISTICS.PACKETS_TO_BS);
subplot(2,2,4);
plot(r,STATISTICS.COUNTCHS);
%STATISTICS���ṹ�����飬���������5��������
%countCHs(r+1��,ÿһ����ѡ���Ĵ�ͷ��Ŀ;
%packets_TO_BS(r+1),��վ�յ������ݰ�����;
%PACKETS_TO_CH(r+1),��ͷ�յ������ݰ�����;
%first_dead,��һ���ڵ�������ʱ��;
%teenth_dead=r,10%�Ľڵ�������ʱ�䣻
%dead(r+1),ÿһ�ֵ������ڵ�����
%allive(r+1),ÿһ�ֵĻ�ڵ�����

