clear;
m=1;
for a=0:0.5:5
%1.��ʼ�����趨ģ��
%.�������ڵ��������(��λ M)
xm=100;
ym=100;
%(1)��۽��������
sink.x=xm;
sink.y=ym;
%�����ڴ�������
n=100;
%��ͷ�Ż���������ѡ��ͷ�ĸ��ʣ�
p=0.1;
P=0.1;
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
rmax=50;
%������� do
do=sqrt(Efs/Emp);
Et=0;
%2.���ߴ���������ģ�Ͳ���ģ��
%�������ߴ���������,�������ھ���Ͷ��100���ڵ�,������ͼ��
for i=1:1:n
    S1(i).xd=rand(1,1)*xm;
    S2(i).xd=S1(i).xd;
    S3(i).xd=S1(i).xd;
    S4(i).xd=S3(i).xd;
    XR4(i)=S4(i).xd;
    XR3(i)=S3(i).xd;
    XR2(i)=S2(i).xd;
    XR1(i)=S1(i).xd;
    S1(i).yd=rand(1,1)*ym;
    S2(i).yd=S1(i).yd;
    S3(i).yd=S1(i).yd;
    S4(i).yd=S3(i).yd;
    YR4(i)=S4(i).yd;
    S4(i).G=0;
    YR3(i)=S3(i).yd;
    S3(i).G=0;
    YR2(i)=S2(i).yd;
    YR1(i)=S1(i).yd;
    S1(i).G=0;
    S2(i).G=0;
    S1(i).E=Eo*(1+rand*a);
    S2(i).E=S1(i).E;
    S3(i).E=S1(i).E;
    S4(i).E=S3(i).E;
    E3(i)= S3(i).E;
    E4(i)= S4(i).E;
    Et=Et+E3(i);
    %initially there are no cluster heads only nodes
    S1(i).type='N';
    S2(i).type='N';
    S3(i).type='N';
    S4(i).type='N';
    
end
S2(n+1).xd=sink.x;
S2(n+1).yd=sink.y;
%3.��������ģ��
%��ͷ�ڵ���
countCHs2=0;
cluster=1;%�˶����Ŀ�Ľ����Ǹ���һ��1��ʼ���±�����������Ĵ�ͷ��Ӧ�û���ȥ1
flag_first_dead2=0;
flag_teenth_dead2=0;
%�����ڵ���
dead2=0;
first_dead2(m)=0;
teenth_dead2(m)=0;
%��ڵ���
allive2=n;
%(1)ѭ��ģʽ�趨
for r=0:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r;
  %ÿ��һ����ת����(������Ϊ10��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣��ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ��
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S2(i).G=0;
        S2(i).cl=0;
    end
  end
%(2)�����ڵ���ģ��
dead2=0;
for i=1:1:n
    %������������ڵ�
    if (S2(i).E<=0)
        dead2=dead2+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead2==1)
           if(flag_first_dead2==0)
              first_dead2(m)=r;
              flag_first_dead2(m)=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead2==0.1*n)
           if(flag_teenth_dead2==0)
              teenth_dead2(m)=r;
              flag_teenth_dead2=1;
           end
        end
    end
    if S2(i).E>0
        S2(i).type='N';
    end
end
%(4)��ͷѡ��ģ��
countCHs2=0;
cluster2=1;
for i=1:1:n
   if(S2(i).E>0)
   temp_rand=rand;     
   if ( (S2(i).G)<=0)  
       %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs2=countCHs2+1;
            S2(i).type='C';
            S2(i).G=round(1/p)-1;
            C2(cluster2).xd=S2(i).xd;
            C2(cluster2).yd=S2(i).yd;
           distance=sqrt( (S2(i).xd-(S2(n+1).xd) )^2 + (S2(i).yd-(S2(n+1).yd) )^2 );
            C2(cluster2).distance=distance;
            C2(cluster2).id=i;
            X2(cluster2)=S2(i).xd;
            Y2(cluster2)=S2(i).yd;
            cluster2=cluster2+1;
           %�����ͷ����4000bit���ݵ���վ���������ģ�����Ӧ�����нڵ������ͷÿһ�ַ���4000bit���ݣ�
           distance;
            if (distance>do)
                S2(i).E=S2(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S2(i).E=S2(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    
    end
    % S2(i).G=S2(i).G-1;  
end 
end
%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
for i=1:1:n
   if ( S2(i).type=='N' && S2(i).E>0 )
     if(cluster2-1>=1)
       min_dis=sqrt( (S2(i).xd-S2(n+1).xd)^2 + (S2(i).yd-S2(n+1).yd)^2 );
       min_dis_cluster=0;
       for c=1:1:cluster2-1
           temp=min(min_dis,sqrt( (S2(i).xd-C2(c).xd)^2 + (S2(i).yd-C2(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       %���ڽڵ㣨����4000bit���ݣ���������
       if(min_dis_cluster~=0)    
            min_dis;
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %��ͷ�����ܺ��ں���һ���ڽڵ�4000bit���ݣ�����������
            S2(C2(min_dis_cluster).id).E = S2(C2(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
        else 
            min_dis;
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        end
        S2(i).min_dis=min_dis;
        S2(i).min_dis_cluster=min_dis_cluster;
    else
        min_dis=sqrt( (S2(i).xd-S2(n+1).xd)^2 + (S2(i).yd-S2(n+1).yd)^2 );
            if (min_dis>do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S2(i).E=S2(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        end
  end
end
end
d1=0.765*xm/2;
K=sqrt(0.5*n*do/pi)*xm/d1^2;
d2=xm/sqrt(2*pi*K);
Er=4000*(2*n*ETX+n*EDA+K*Emp*d1^4+n*Efs*d2^2);
S4(n+1).xd=sink.x;
S4(n+1).yd=sink.y;
countCHs4=0;
cluster=1;%�˶����Ŀ�Ľ����Ǹ���һ��1��ʼ���±�����������Ĵ�ͷ��Ӧ�û���ȥ1
flag_first_dead4=0;
flag_teenth_dead4=0;
%�����ڵ���
dead4=0;
first_dead4(m)=0;
teenth_dead4(m)=0;
%��ڵ���
allive4=n;
%(1)ѭ��ģʽ�趨
for r=0:1:rmax     %�� for ѭ������������г���������ڣ�ֱ�����һ end �Ž���ѭ��
    r;
  %ÿ��һ����ת����(������Ϊ10��)ʹ���ڵ��S��i��.G�������ò������ں���Ĵ�ѡ�٣��ڸ���ת�������ѵ�ѡ����ͷ�Ľڵ㲻���ٵ�ѡ���ָ�Ϊ��
  if(mod(r, round(1/P) )==0)
    for i=1:1:n
        S4(i).G=0;
        S4(i).cl=0;
    end
  end
Ea=Et*(1-r/rmax)/n;
%(2)�����ڵ���ģ��
dead4=0;
for i=1:1:n
    %������������ڵ�
    if (S4(i).E<=0)
        dead4=dead4+1; 
        %(3)��һ�������ڵ�Ĳ���ʱ��(���ִα�ʾ)
        %��һ���ڵ�����ʱ��
        if (dead4==1)
           if(flag_first_dead4==0)
              first_dead4(m)=r;
              flag_first_dead4=1;
           end
        end
        %10%�Ľڵ�����ʱ��
        if(dead4==0.1*n)
           if(flag_teenth_dead4==0)
              teenth_dead4(m)=r;
              flag_teenth_dead4=1;
           end
        end

    end
    if S4(i).E>0
        S4(i).type='N';
    end
end
%(4)��ͷѡ��ģ��
countCHs4=0;
cluster4=1;
for i=1:1:n
 if Ea>0
 p(i)=P*n*S4(i).E*E4(i)/(Et*Ea);
 if(S4(i).E>0)
   temp_rand=rand;     
   if ( (S4(i).G)<=0)  
       %��ͷ��ѡ�٣���ѡ�Ĵ�ͷ��Ѹ�������Ŵ�����������������ı�����
        if(temp_rand<= (p(i)/(1-p(i)*mod(r,round(1/p(i))))))
            countCHs4=countCHs4+1;
             S4(i).type='C';
            S4(i).G=round(1/p(i))-1;
            C4(cluster4).xd=S4(i).xd;
            C4(cluster4).yd=S4(i).yd;
           distance=sqrt( (S4(i).xd-(S4(n+1).xd) )^2 + (S4(i).yd-(S4(n+1).yd) )^2 );
            C4(cluster4).distance=distance;
            C4(cluster4).id=i;
            X4(cluster4)=S4(i).xd;
            Y4(cluster4)=S4(i).yd;
            cluster4=cluster4+1;
           %�����ͷ����4000bit���ݵ���վ���������ģ�����Ӧ�����нڵ������ͷÿһ�ַ���4000bit���ݣ�
           distance;
            if (distance>do)
                S4(i).E=S4(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance )); 
            end
            if (distance<=do)
                S4(i).E=S4(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance )); 
            end
        end     
    end
    % S4(i).G=S4(i).G-1;  
end 
 end
end
%(5)���ڳ�Աѡ���ͷģ��(���ص��γ�ģ��)
%���ڳ�Ա�Դ�ͷ��ѡ�񣨼��ص��γɣ��㷨
for i=1:1:n
   if ( S4(i).type=='N' && S4(i).E>0 )
     if(cluster4-1>=1)
       min_dis=sqrt( (S4(i).xd-S4(n+1).xd)^2 + (S4(i).yd-S4(n+1).yd)^2 );
       min_dis_cluster=0;
       for c=1:1:cluster4-1
           temp=min(min_dis,sqrt( (S4(i).xd-C4(c).xd)^2 + (S4(i).yd-C4(c).yd)^2 ) );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       %���ڽڵ㣨����4000bit���ݣ���������
       if(min_dis_cluster~=0)    
            min_dis;
            if (min_dis>do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        %��ͷ�����ܺ��ں���һ���ڽڵ�4000bit���ݣ�����������
            S4(C4(min_dis_cluster).id).E = S4(C4(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 
        else 
            min_dis;
            if (min_dis>do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        end
        S4(i).min_dis=min_dis;
       S4(i).min_dis_cluster=min_dis_cluster;
   else
            min_dis=sqrt( (S4(i).xd-S4(n+1).xd)^2 + (S4(i).yd-S4(n+1).yd)^2 );
            if (min_dis>do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S4(i).E=S4(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
        end
  end
end
end
m=m+1;
end
r=0:5:50;
subplot(2,1,1);
plot(r,first_dead2,'-r',r,first_dead4,'-.m');
legend('leachm2','DEEC2');
xlabel('x(time)');
ylabel('y(dead)');
title('\bf leachm2��DEEC2)���������ڶԱ�');
subplot(2,1,2);
plot(r,teenth_dead2,'-r',r,teenth_dead4,'-.m');
legend('leachm2','DEEC2');
xlabel('x(time)');
ylabel('y(data)');
title('\bf leachm2��DEEC2���ȶ����ڶԱ�');