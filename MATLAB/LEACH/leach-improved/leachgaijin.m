%Field Dimensions - x and y maximum (in meters)
xm=300;
ym=300;
%x and y Coordinates of the Sink
sink.x=0.5*xm;
sink.y=0.5*ym;
%Number of Nodes in the field
n=900;
%Optimal Election Probability of a node
%to become cluster head
p=0.1;
%Energy Model (all values in Joules)
%Initial Energy 
Eo=0.5;
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
m=0.1;
%\alpha
a=1;
cc=10;
CM=32;
DM=4000;
%maximum number of rounds
rmax=100;
%%%%%%%%%%%%%%%%%%%%%%%%% END OF PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%Computation of do
do=sqrt(Efs/Emp);
%Creation of the random Sensor Network
figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    %initially there are no cluster heads only nodes
    S(i).type='N';
   
    temp_rnd0=i;
    %Random Election of Normal Nodes
    if (temp_rnd0>=m*n+1) 
        S(i).E=Eo;
        S(i).ENERGY=0;
       % plot(S(i).xd,S(i).yd,'o');
       % hold on;
    end
    %Random Election of Advanced Nodes
    if (temp_rnd0<m*n+1)  
        S(i).E=Eo*(1+a);
        S(i).ENERGY=1;
       % plot(S(i).xd,S(i).yd,'+');
       % hold on;
    end
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%plot(S(n+1).xd,S(n+1).yd,'x');
    
        
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
for r=0:1:rmax
    r+1
  %Operation for epoch
  if(mod(r, round(1/p) )==0)
    for i=1:1:n
        S(i).G=0;
        %S(i).cl=0;
    end
  end
hold off;
EJ(r+1)=0;
%Number of dead nodes
dead=0;
%Number of dead Advanced Nodes
dead_a=0;
%Number of dead Normal Nodes
dead_n=0;
%counter for bit transmitted to Bases Station and to Cluster Heads
%packets_TO_BS=0;
%packets_TO_CH=0;
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
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    else
        EJ(r+1)=EJ(r+1)+S(i).E;
        S(i).type='N';
     %   if (S(i).ENERGY==0)  
     %   plot(S(i).xd,S(i).yd,'o');
      %  end
      %  if (S(i).ENERGY==1)  
       % plot(S(i).xd,S(i).yd,'+');
     %   end
      %  hold on;
    end
end
%plot(S(n+1).xd,S(n+1).yd,'x');
%nl=n-dead;
%p=do*sqrt(xm*ym/(2*pi*nl))/(2/3*xm*ym-(sink.x+sink.y)*sqrt(xm*ym)+sink.x^2+sink.y^2)
p=1/(cc+1);

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;

%When the first node dies
if (dead==1)
    if(flag_first_dead==0)
        first_dead=r
        flag_first_dead=1;
    end
end
countCHs=0;
cluster=1;
for i=1:1:n
   if(S(i).E>0)
   temp_rand=rand;     
   if ( (S(i).G)<=0)
%Election of Cluster Heads
         if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
            countCHs=countCHs+1;
                        
            S(i).type='C';
            S(i).G=round(1/p)-1;
            C(cluster).xd=S(i).xd;
            C(cluster).yd=S(i).yd;
            %plot(S(i).xd,S(i).yd,'k*');
            
            distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
            C(cluster).distance=distance;
            C(cluster).id=i;
            X(cluster)=S(i).xd;
            Y(cluster)=S(i).yd;
            packets_TO_BS(cluster)=1;
            cluster=cluster+1;
            
            %Calculation of Energy dissipated
            
        end     
    
    end
  end 
end
STATISTICS(r+1).CLUSTERHEADS=cluster-1;
CLUSTERHS(r+1)=cluster-1;
%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 )
     if(cluster-1>=1)
       min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
       min_dis_cluster=1;
       for c=1:1:cluster-1
           temp=sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 );
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
       end
       packets_TO_BS(min_dis_cluster)=packets_TO_BS(min_dis_cluster)+1;
       %line([S(i).xd, C(min_dis_cluster).xd],[S(i).yd, C(min_dis_cluster).yd]);
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            Er1=ERX*CM*(cluster+1);
            if (min_dis>do)
                Et1=ETX*(CM+DM)+Emp*(CM+DM)* min_dis * min_dis * min_dis * min_dis; 
            end
            if (min_dis<=do)
                Et1=ETX*(CM+DM)+Efs*(CM+DM)*min_dis * min_dis; 
            end
            S(i).E=S(i).E-Er1-Et1;
            EJ(r+1)=EJ(r+1)-Er1-Et1;
       S(i).min_dis=min_dis;
       S(i).min_dis_cluster=min_dis_cluster;
     end   
   end
end
%a=zeros(cluster-1,1)
for c=1:1:cluster-1
    CEr1=ERX*(CM+DM)*(packets_TO_BS(c)-1);
    distemp=0;
    for pinkx=0:xm:xm
        for pinky=0:ym:ym
         dispink=sqrt( (X(c)-pinkx)^2 + (Y(c)-pinky)^2 );   
         if(dispink>distemp) 
             distemp=dispink;
         end
        end
    end
   distbroadcast(c)= distemp;
   if (distbroadcast(c)>do)
         CEt1=ETX*CM+Emp*CM*distbroadcast(c)*distbroadcast(c)*distbroadcast(c)*distbroadcast(c);
   end
   if (distbroadcast(c)<=do)
         CEt1=ETX*CM+Efs*CM*distbroadcast(c)*distbroadcast(c);
   end
  
   S(C(c).id).E=S(C(c).id).E-CEr1-CEt1;
   
   if (S(C(c).id).E<=0)
       % plot(S(i).xd,S(i).yd,'red .');
        dead=dead+1;
        if(S(i).ENERGY==1)
            dead_a=dead_a+1;
        end
        if(S(i).ENERGY==0)
            dead_n=dead_n+1;
        end
        hold on;    
    else
         EJ(r+1)=EJ(r+1)-CEr1-CEt1;
     
    end
end

STATISTICS(r+1).DEAD=dead;
DEAD(r+1)=dead;
DEAD_N(r+1)=dead_n;
DEAD_A(r+1)=dead_a;


   Et=0;
   Et=Et+S(C(c).id).E;
  % EJ(r+1)=EJ(r+1)-CEr1-CEt1;
  %PACKETS_TO_BS(r+1)=packets_TO_BS(c);
hold on;

%E=EJ(r+1)/nl;
%b=zeros(,1)
%for i=1:1:cluster-1
%   if(a(i)>E)
  %      b(i)=a(i);
   % en
%end
CH=1;
while(cluster-1>20)%contiue make cluster
       %select sec cluster heads
       E=Et/(cluster-1);
for c=1:1:cluster-1
    if (S(C(c).id).E>E)
       C(c).type='D';
       D(CH).xd=X(c);
       D(CH).yd=Y(c);
       distance=sqrt( (X(c)-S(n+1).xd)^2 + (Y(c).yd-S(n+1).yd)^2 );
       D(CH).distance=distance;
       D(CH).id=C(c).id;
       X(CH)= C(cluster).xd;
       Y(CH)= C(cluster).yd;
       packets_TO_BS(CH)=1;
       CH=CH+1; 
    end
end
       
STATISTICS(r+1).CLUSTERHEADS=CH-1;
CLUSTERHS(r+1)=CH-1;

for c=1:1:cluster-1
   if (  C(c).type=='C' && S(C(c).id).E>0 )
     if(CH-1>=1)
      % min_dis=sqrt( (C(c).xd-S(n+1).xd)^2 + (C(c).yd-S(n+1).yd)^2 );
       min_dis_CH=1;
       for d=1:1:CH-1
           temp=sqrt( (D(d).xd-C(c).xd)^2 + (D(d).yd-C(c).yd)^2 );
            if ( temp<C(c).distance )
               min_dis=temp;
               min_dis_CH=d;
           end
       end
       packets_TO_BS(min_dis_CH)=packets_TO_BS(min_dis_CH)+1;
       %line([S(i).xd, C(min_dis_cluster).xd],[S(i).yd, C(min_dis_cluster).yd]);
       
       %Energy dissipated by associated Cluster Head
            min_dis;
            Er1=ERX*CM*(CH+1);
            if (min_dis>do)
                Et1=ETX*(CM+DM)+Emp*(CM+DM)* min_dis * min_dis * min_dis * min_dis; 
            end
            if (min_dis<=do)
                Et1=ETX*(CM+DM)+Efs*(CM+DM)*min_dis * min_dis; 
            end
            S(C(c).id).E=S(C(c).id).E-Er1-Et1;
            EJ(r+1)=EJ(r+1)-Er1-Et1;
       S(C(c).id).min_dis=min_dis;
       S(C(c).id).min_dis_CH=min_dis_CH;
     end   
   end
end

for d=1:1:CH-1
    CEr1=ERX*(CM+DM)*(packets_TO_BS(d)-1);
   distbroadcast(d)= D(d).distance;
   if (distbroadcast(d)>do)
         CEt1=ETX*CM+Emp*CM*distbroadcast(d)*distbroadcast(d)*distbroadcast(d)*distbroadcast(d);
   end
   if (distbroadcast(d)<=do)
         CEt1=ETX*CM+Efs*CM*distbroadcast(d)*distbroadcast(d);
   end
  
   if(packets_TO_BS(d)<=cc)
         l=1;
   else
       l=ceil(packets_TO_BS(d)/cc);
   end
   
   if (D(d).distance>do)
          CEt2=(ETX+EDA)*DM*l+ Emp*DM*l* D(d).distance*D(d).distance*D(d).distance*D(d).distance ; 
   end
   if (D(d).distance<=do)
          CEt2=(ETX+EDA)*DM*l+ Efs*DM*l* D(d).distance*D(d).distance ; 
   end
   
   S(D(d).id).E=S(D(d).id).E-CEr1-CEt1-CEt2;
   EJ(r+1)=EJ(r+1)-CEr1-CEt1-CEt2;
  PACKETS_TO_BS(r+1)=packets_TO_BS(d);
hold on;
end

end
end
%subplot(2,2,1);
%r=0:1:rmax;
%plot(r,DEAD(r+1));
%box off;
%title('图一：死亡节点数统计图');
%xlabel('工作轮数（/轮）');
%ylabel('死亡节点数（/个）');
%subplot(2,2,2);
r=0:1:rmax;
plot(r,PACKETS_TO_BS(r+1));
box off;
title('图二：基站接收数据包统计图');
xlabel('工作轮数（/轮）');
ylabel('数据包个数（/个）');
%subplot(2,2,3);
%r=0:1:rmax;
%plot(r,EJ(r+1));
%box off;
%title('图三：网络剩余能量统计图');
%xlabel('工作轮数（/轮）');
%ylabel('网络剩余能量（/焦耳）');
