clear;
%1.��ʼ�����趨ģ��
%.�������ڵ��������(��λ M)
xm=30;
ym=30;
%(1)��۽��������
sink.x=xm;
sink.y=ym;
%�����ڴ�������
n=30;
%��ͷ�Ż���������ѡ��ͷ�ĸ��ʣ�
p=0.3;
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
%���ѭ������
rmax=1;
%������� do
do=sqrt(Efs/Emp);
%2.���ߴ���������ģ�Ͳ���ģ��
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
%3.��������ģ��
%��ͷ�ڵ���
countCHs1=0;
cluster1=1;%�˶����Ŀ�Ľ����Ǹ���һ��1��ʼ���±�����������Ĵ�ͷ��Ӧ�û���ȥ1
flag_first_dead1=0;
flag_all_dead1=0;
%�����ڵ���
dead1=0;
first_dead1=0;
all_dead1=0;
%��ڵ���
allive1=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS1=0;
packets_TO_CH1=0;
%(1)ѭ��ģʽ�趨
for r=0:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��    
  %ÿ��һ����ת����ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣��ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ��
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
        dead1=dead1+1; 
        plot(S(i).xd,S(i).yd,'g+');
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
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
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            plot(S(i).xd,S(i).yd,'k*');
            C(cluster).id=i;
            X(cluster)=S(i).xd;%CH��x����
            Y(cluster)=S(i).yd;%CH��y����
            cluster=cluster+1;
        end
   end
   end
end
STATISTICS.COUNTCHS(r+1)=countCHs;
P=[X S(n+1).xd];
Q=[Y S(n+1).yd];
%����do����Ȩ�ؾ���
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


%��·��ת��
[d,index1,index2]=short(D);

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



