clear;
xm=200;
ym=200;
Eelec = 50*10^(-9); % energy consumption per bit
%Eamp = 100*10^(-12); % multipath fading
%EDA = 5*10^(-9); % energy for data aggregation
n=100;
Eo=0.5;          %energy supplied to each node
ch=n/10;
ETX=50*0.000000001;     %transmiter energy per node
ERX=50*0.000000001;        %reciever energy per mode

Efs=10*0.000000000001;     %amplification energy when d is less than d0
Emp=0.0013*0.000000000001;      %amplification energy  when d is greater than d0
Efs1=Efs/10;   % amp energy just for intra cluster communication.
Emp1=Emp/10;
%Data Aggregation Energy
EDA=5*0.000000001;

a=Eo/2;                %?
%temprature range
tempi=50;
tempf=200;

%Thresholod for transmiting data to the cluster head
h=100;                                 %%%%%%Hard Thres%%%%hold H(t)
s=2;                                   %%%%%%Soft thres%%%%hold  S(t)

sv=0;                                  %%%%%%previously Sensed value S(v)


do=sqrt(Efs/Emp);       %distance between cluster head and base station
do1=sqrt(Efs1/Emp1);    

%x=0;
%cluster=0;
bs.xd=xm*0.5;
bs.yd=ym*0.5;
rmax= 100;
Eo=0.5;  %energy of each node
p=0.1; %probability of node to become cluster heads
k=4096; %packet length=4Bytes
figure(1); %random node distribution in field
for i=1:1:n
    S(i).xd=rand(1,1)*xm; 
    S1(i).xd1=rand(1,1)*xm;
    S(i).yd=rand(1,1)*ym;
    S1(i).yd1=rand(1,1)*ym;
    
    S(i).G=0;
    S1(i).G=0;
    S(i).type='N';
    S1(i).type1='N';
    %Normal Nodes 
    S(i).E=Eo;
    S1(i).E1=Eo;
    
    plot(S(i).xd,S(i).yd,'o');
    %plot(S1(i).xd1,S1(i).yd1,'o');
   
    hold on;
end
S(n+1).xd=bs.xd;
S(n+1).yd=bs.yd; %base station(sink)
S1(n+1).xd1=bs.xd; %sink is a n+1 node, x-axis postion of a node
S1(n+1).yd1=bs.yd; %sink is a n+1 node, y-axis postion of a node

plot(S(n+1).xd,S(n+1).yd,'g+');
%plot(S1(n+1).xd1,S1(n+1).yd1,'+r');

cluster=1;
x=0;
allive=n;
for r=0:1:rmax
    r;
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S(i).G=0;%tells us if node has become cluster head in last 1/p rounds
        end
    end

    hold off;


%Number of dead nodes
    dead=0;
   
    %figure(2);

    for i=1:1:n
        if (S(i).E<=0)

%checking if there is a dead node or active node by comparing the energy of each node

            plot(S(i).xd,S(i).yd,'red .');
            dead=dead+1;
            hold on;    
        end
        if(S(i).E>0)
            S(i).type='N'; 
            plot(S(i).xd,S(i).yd,'mo');
        end
            hold on;
    end

 plot(S(n+1).xd,S(n+1).yd,'gx');
hold on;
 STATISTICS.DEAD(r+1)=dead;
 STATISTICS.ALLIVE(r+1)=allive-dead;

%DEAD(r+1)=dead;
%ALLIVE(r+1)=allive-dead;
%When the first node dies
    if (dead==1)
        if(x==0)
        first_dead=r;
        x=1;
        end
    end
 
%count Cluster Heads=0;
cluster=1;
    for i=1:1:n 
        if(S(i).E>0)
        temp_rand=rand;     
            if ( (S(i).G<=0))
 %Election of Cluster Heads using the treshold equation
                if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                S(i).type='C';
                S(i).G=round(1/p)-1;
                C(cluster).xd=S(i).xd;
                C(cluster).yd=S(i).yd;
                plot(S(i).xd,S(i).yd,'*');
                distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
                C(cluster).distance=distance;
                C(cluster).id=i;
                cluster=cluster+1;
                Etx(i)=(Eelec+Emp*distance*distance)*k;
                S(i).E=S(i).E-Etx(i)-EDA; 
                end     
            end
        end 
    end
    
    for i=1:1:n
        if ( S(i).type=='N' && S(i).E>0 )
            if(cluster-1>=1)
                min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                min_dis_cluster=1;
                for c=1:1:cluster-1
                temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                    if ( temp<min_dis )
                        min_dis=temp;
                        min_dis_cluster=c;
                    end
                end

 %calucating energies

                Etx(i)=(Eelec+Emp*min_dis*min_dis)*k;
                S(i).E=S(i).E-Etx(i);  
                if(min_dis>0)
                Erx=Eelec*k;
                S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- (Erx + EDA); 
                end
                S(i).min_dis=min_dis;
                S(i).min_dis_cluster=min_dis_cluster;
            end
        end
    end
hold on;
energy=0;
for i=1:n
    energy=energy+S(i).E;
end
residue_energy(r+1)=energy;
end

for i=1:1:rmax+1
    q(i)=i;
end
STATISTICS.DEAD(r+1);
STATISTICS.ALLIVE(r+1);


% mod leach 

countCHs1=0;         %the number of Stateflow objects in the current context.
cluster1=1;              %first cluster is selected
flag_first_dead1=0;         
flag_teenth_dead1=0;
flag_all_dead1=0;

dead1=0;
first_dead=0;
teenth_dead=0;
all_dead=0;

allive1=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS1=0;
packets_TO_CH1=0;
for r=0:1:rmax
    r;
   cv = tempi + (tempf-tempi).*rand(1,1);  %%%%%%Current sensing value C(v)
   if(mod(r, round(1/p) )==0) %remainder
   for i=1:1:n
       S1(i).G1=0;            % it will assign to the nodes that have not been cluster head .
       %%S(i).cl=0;
   end
   end
%hold off;
dead1=0;

%%Checking dead nodes
% figure(1);
for i=1:1:n

   if (S1(i).E1<=0)
       dead1=dead1+1;
       %%plot(S1(i).xd1,S1(i).yd1,'red.')   
       if (dead1==1)
          if(flag_first_dead1==0)
             first_dead1=r;
             flag_first_dead1=1;
          end
       end

       if(dead1==0.1*n)
          if(flag_teenth_dead1==0)
             teenth_dead1=r;
             flag_teenth_dead1=1;
          end
       end
       if(dead1==n)
          if(flag_all_dead1==0)
             all_dead1=r;
             flag_all_dead1=1;
          end
       end
   end
   if S1(i).E1>0
       S1(i).type1='N';
       %plot(S1(i).xd1,S1(i).yd1,'o');
   end
   %hold on;
end
STATISTICS1.DEAD1(r+1)=dead1;
STATISTICS1.ALLIVE1(r+1)=allive1-dead1;

countCHs1=0;
cluster1=1;

if   S1(i).type1=='C' && S1(i).E1>a
for j=1:1:ch
           countCHs1=countCHs1+1;
           S1(i).type1='C';
           %plot(S1(i).xd1,S1(i).yd1,'r*');
           %hold on;
           S1(i).G1=round(1/p)-1;
           C1(cluster1).xd1=S1(i).xd1;
           C1(cluster1).yd1=S1(i).yd1;
           distance1=sqrt( (S1(i).xd1-(S1(n+1).xd1) )^2 + (S1(i).yd1-(S1(n+1).yd1) )^2 );
           C1(cluster1).distance1=distance1;
           C1(cluster1).id=i;
           X1(cluster1)=S1(i).xd1;
           Y1(cluster1)=S1(i).yd1;
           cluster1=cluster1+1;
           distance1;
           
    %   if (cv >= h)
 %test = cv-sv;
 %if (test >= s)
           if (distance1>do)
               S1(i).E1=S1(i).E1- ( (ETX+EDA)*(4000) + Emp*4000*(distance1*distance1*distance1*distance1 ));
           end
           if (distance1<=do)
               S1(i).E1=S1(i).E1- ( (ETX+EDA)*(4000)  + Efs*4000*(distance1 * distance1 ));
           end
           %end
           
           %packets_TO_BS=packets_TO_BS+1;
           %PACKETS_TO_BS(r+1)=packets_TO_BS;
            %          packets_TO_CH=packets_TO_CH+1;
end
else
for i=1:1:n        
  if(S1(i).E1>0)
  temp_rand1=rand;
  if ( (S1(i).G1)<=0)

       if(temp_rand1<= (p/(1-p*mod(r,round(1/p)))))
           countCHs1=countCHs1+1;
           %plot(S1(i).xd1,S1(i).yd1,'r*'); 
           %hold on;
           packets_TO_BS1=packets_TO_BS1+1;
           PACKETS_TO_BS1(r+1)=packets_TO_BS1;
           S1(i).type1='C';
           S1(i).G1=round(1/p)-1;
           C1(cluster1).xd1=S1(i).xd1;
           C1(cluster1).yd1=S1(i).yd1;
           distance1=sqrt( (S1(i).xd1-(S1(n+1).xd1) )^2 + (S1(i).yd1-(S1(n+1).yd1) )^2 );
           C1(cluster1).distance1=distance1;
           C1(cluster1).id1=i;
           X1(cluster1)=S1(i).xd1;
           Y1(cluster1)=S1(i).yd1;
           cluster1=cluster1+1;
             
    %   if (cv >= h)
 %test = cv-sv;
 %if (test >= s)
               
           distance1;
           if (distance1>do)
               S1(i).E1=S1(i).E1- ( (ETX+EDA)*(4000) + Emp*4000*(distance1*distance1*distance1*distance1 ));
           end
           if (distance1<=do)
               S1(i).E1=S1(i).E1- ( (ETX+EDA)*(4000)  + Efs*4000*(distance1 * distance1 ));
           end
       end
  end
  % end
   % S(i).G=S(i).G-1;

  end
end
end
STATISTICS1.COUNTCHS1(r+1)=countCHs1;

for i=1:1:n
  if ( S1(i).type1=='N' && S1(i).E1>0 )
    if(cluster1-1>=1)
      min_dis1=Inf;
      min_dis_cluster1=0;
      for c1=1:1:cluster1-1
          temp1=min(min_dis1,sqrt( (S1(i).xd1-C1(c1).xd1)^2 + (S1(i).yd1-C1(c1).yd1)^2 ) );
          if ( temp1<min_dis1 )
              min_dis1=temp1;
              min_dis_cluster1=c1;
          end
      end
% if (cv >= h)
 %test = cv-sv;
 %if (test >= s)
           min_dis1;
           if (min_dis1>do1)
               S1(i).E1=S1(i).E1- ( ETX*(4000) + Emp1*4000*( min_dis1 *min_dis1 * min_dis1 * min_dis1));
           end
          if (min_dis1<=do1)
               S1(i).E1=S1(i).E1- ( ETX*(4000) + Efs1*4000*( min_dis1 * min_dis1));
           end
                  S1(C1(min_dis_cluster1).id1).E1 =S1(C1(min_dis_cluster1).id1).E1- ( (ERX + EDA)*4000 );
           packets_TO_CH1=packets_TO_CH1+1;
 %end
  %sv
        
       S1(i).min_dis1=min_dis1;
       S1(i).min_dis_cluster1=min_dis_cluster1;
   else
       min_dis1=sqrt( (S1(i).xd1-S1(n+1).xd1)^2 + (S1(i).yd1-S1(n+1).yd1)^2 );
           if (min_dis1>do)
               S1(i).E1=S1(i).E1- ( ETX*(4000) + Emp*4000*( min_dis1 *min_dis1 * min_dis1 * min_dis1));
           end
           if (min_dis1<=do)
               S1(i).E1=S1(i).E1- ( ETX*(4000) + Efs*4000*( min_dis1 * min_dis1));
           end
           packets_TO_BS1=packets_TO_BS1+1;
    
  sv=cv;
    end
 end
end
%STATISTICS1.PACKETS_TO_CH1(r+1)=packets_TO_CH1;
%STATISTICS1.PACKETS_TO_BS1(r+1)=packets_TO_BS1;
end
STATISTICS1.DEAD1(r+1);
STATISTICS1.ALLIVE1(r+1);

r=0:rmax;
figure(2);
plot(r,STATISTICS.DEAD(r+1),'-r',r,STATISTICS1.DEAD1(r+1),'-g')
grid on
legend('LEACH','LEACH-C');
xlabel('Rounds');
ylabel('Dead Nodes');
title('Dead Nodes');
figure(3);
plot(r,STATISTICS.ALLIVE(r+1),'-r',r,STATISTICS1.ALLIVE1(r+1),'-g')
grid on
legend('LEACH','LEACH-C');
xlabel('Rounds');
ylabel('Allive Nodes');
title('Allive nodes');

%figure(3);
%plot(q,DEAD);
%title('dead nodes per round')
%figure(4);
%plot(q,residue_energy);

%title('residual energy per round')
%figure(5);
%plot(q,alive);

%title('nodes alive per round')
%DEAD(rmax)

