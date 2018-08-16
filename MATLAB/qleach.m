clear;

xm=100;
ym=100;
sink.x=0.5*xm;
sink.y=0.5*ym;
n=100;
p=0.05;
Eo=0.5;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
EDA=5*0.000000001;
rmax=2500;
do=sqrt(Efs/Emp);

figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    S(i).type='N';
    S(i).E=Eo;
    plot(S(i).xd,S(i).yd,'ko');
    hold on;
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'rx');

S5=S;

% Q-LEACH,4个分区，4个CH点
S1=[];
S2=[];
S3=[];
S4=[];

for i=1:1:n
    if S(i).xd<=50
        if S(i).yd<=50          
          S1=[S1,i];
        else
          S2=[S2,i];
        end
    else
        if S(i).yd<=50
            S3=[S3,i];
        else
            S4=[S4,i];
        end
    end
end

figure(2);
for i=1:1:n    
        if ismember(i,S1)            
            plot(S(i).xd,S(i).yd,'ko');
            hold on;
        elseif ismember(i,S2)
            plot(S(i).xd,S(i).yd,'ro');
            hold on;
        elseif ismember(i,S3)
            plot(S(i).xd,S(i).yd,'bo');
            hold on;
        else 
            plot(S(i).xd,S(i).yd,'yo');
            hold on;
        end     
end

S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'rx');


% figure(2);

for r=0:rmax
    r
    
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S(i).G=0;            
        end
    end
    
    
    dead=0;
    
    packets_TO_BS=0;
    packets_TO_CH=0;
    PACKETS_TO_CH(r+1)=0;
    PACKETS_TO_BS(r+1)=0;
    
    for i=1:1:n
        if (S(i).E<=0)
            dead=dead+1;
        else
            S(i).type='N';
        end       
    end
    
    STATISTICS(r+1).DEAD=dead;
    DEAD(r+1)=dead;
    
%     if (dead==1)
%         if(flag_first_dead==0)
%             first_dead=r;
%             flag_first_dead=1;
%         end
%     end    
    
countCHs1=0;
cluster1=1;
countCHs2=0;
cluster2=1;
countCHs3=0;
cluster3=1;
countCHs4=0;
cluster4=1;
    
    for i=1:1:n
        if ismember(i,S1)             
            if countCHs1<=4
                if(S(i).E>0)
                    temp_rand=rand;
                    if ( (S(i).G)<=0)
                        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                            countCHs1=countCHs1+1;
%                             packets_TO_BS1=packets_TO_BS1+1;
%                             PACKETS_TO_BS1(r+1)=packets_TO_BS1;
                            
                            S(i).type='C';
                            S(i).G=round(1/p)-1;
                            C1(cluster1).xd=S(i).xd; %存S1中的CH坐标
                            C1(cluster1).yd=S(i).yd;
                            
                            distance1=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%TO BS的距离
                            C1(cluster1).distance=distance1;
                            C1(cluster1).id=i;
                            X1(cluster1)=S(i).xd;
                            Y1(cluster1)=S(i).yd;
                            cluster1=cluster1+1;
                            
                            %Calculation of Energy dissipated
                            if (distance1>do)
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance1*distance1*distance1*distance1 ));
                            else
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance1 * distance1 ));
                            end
                            
                            STATISTICS(r+1).CLUSTERHEADS1=cluster1-1;
                            CLUSTERHS1(r+1)=cluster1-1;
                        end
                    end
                end
            end
        elseif ismember(i,S2) 
            if countCHs2<=4
                if(S(i).E>0)
                    temp_rand=rand;
                    if ( (S(i).G)<=0)
                        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                            countCHs2=countCHs2+1;
%                             packets_TO_BS2=packets_TO_BS2+1;
%                             PACKETS_TO_BS2(r+1)=packets_TO_BS2;
                            
                            S(i).type='C';
                            S(i).G=round(1/p)-1;
                            C2(cluster2).xd=S(i).xd; %存S2中的CH坐标
                            C2(cluster2).yd=S(i).yd;
                            
                            distance2=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%TO BS的距离
                            C2(cluster2).distance=distance2;
                            C2(cluster2).id=i;
                            X2(cluster2)=S(i).xd;
                            Y2(cluster2)=S(i).yd;
                            cluster2=cluster2+1;
                            
                            %Calculation of Energy dissipated
                            if (distance2>do)
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance2*distance2*distance2*distance2 ));
                            else
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance2 * distance2 ));
                            end
                            
                            STATISTICS(r+1).CLUSTERHEADS2=cluster2-1;
                            CLUSTERHS2(r+1)=cluster2-1;
                        end
                    end
                end
            end    
            
            
        elseif ismember(i,S3) 
            if countCHs3<=4
                if(S(i).E>0)
                    temp_rand=rand;
                    if ( (S(i).G)<=0)
                        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                            countCHs3=countCHs3+1;
%                             packets_TO_BS3=packets_TO_BS3+1;
%                             PACKETS_TO_BS3(r+1)=packets_TO_BS3;
                            
                            S(i).type='C';
                            S(i).G=round(1/p)-1;
                            C3(cluster3).xd=S(i).xd; %存S3中的CH坐标
                            C3(cluster3).yd=S(i).yd;
                            
                            distance3=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%TO BS的距离
                            C3(cluster3).distance=distance3;
                            C3(cluster3).id=i;
                            X3(cluster3)=S(i).xd;
                            Y3(cluster3)=S(i).yd;
                            cluster3=cluster3+1;
                            
                            %Calculation of Energy dissipated
                            if (distance3>do)
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance3*distance3*distance3*distance3 ));
                            else
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance3 * distance3 ));
                            end
                            
                            STATISTICS(r+1).CLUSTERHEADS3=cluster3-1;
                            CLUSTERHS3(r+1)=cluster3-1;
                        end
                    end
                end
            end
        else
            if countCHs4<=4
                if(S(i).E>0)
                    temp_rand=rand;
                    if ( (S(i).G)<=0)
                        if(temp_rand<= (p/(1-p*mod(r,round(1/p)))))
                            countCHs4=countCHs4+1;
%                             packets_TO_BS4=packets_TO_BS4+1;
%                             PACKETS_TO_BS4(r+1)=packets_TO_BS4;
                            
                            S(i).type='C';
                            S(i).G=round(1/p)-1;
                            C4(cluster4).xd=S(i).xd; %存S1中的CH坐标
                            C4(cluster4).yd=S(i).yd;
                            
                            distance4=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%TO BS的距离
                            C4(cluster4).distance=distance4;
                            C4(cluster4).id=i;
                            X4(cluster4)=S(i).xd;
                            Y4(cluster4)=S(i).yd;
                            cluster4=cluster4+1;
                            
                            %Calculation of Energy dissipated
                            if (distance4>do)
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance4*distance4*distance4*distance4 ));
                            else
                                S(i).E=S(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance4 * distance4 ));
                            end
                            
                            STATISTICS(r+1).CLUSTERHEADS4=cluster4-1;
                            CLUSTERHS4(r+1)=cluster4-1;
                        end
                    end
                end
            end
        end
    end
    
    for i=1:1:n        
        if ( S(i).type=='N' && S(i).E>0 )
            if ismember(i,S1)
                if(cluster1-1>=1)
                    min_dis1=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                    min_dis_cluster1=1;
                    for c1=1:1:cluster1-1
                        temp1=min(min_dis1,sqrt( (S(i).xd-C1(c1).xd)^2 + (S(i).yd-C1(c1).yd)^2 ) );
                        if ( temp1<min_dis1 )
                            min_dis1=temp1;
                            min_dis_cluster1=c1;
                        end
                    end
                    if (min_dis1>do)
                        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis1 * min_dis1 * min_dis1 * min_dis1));
                    else
                        S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis1 * min_dis1));
                    end
                    if(min_dis1>0)
                        S(C1(min_dis_cluster1).id).E = S(C1(min_dis_cluster1).id).E- ( (ERX + EDA)*4000 ); 
%                         PACKETS_TO_CH1(r+1)=n-dead-cluster1+1;
                    end
                    S(i).min_dis1=min_dis1;
                    S(i).min_dis_cluster1=min_dis_cluster1;
                end
                
            elseif ismember(i,S2)
                if(cluster2-1>=1)
                    min_dis2=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                    min_dis_cluster2=1;
                    for c2=1:1:cluster2-1
                        temp2=min(min_dis2,sqrt( (S(i).xd-C2(c2).xd)^2 + (S(i).yd-C2(c2).yd)^2 ) );
                        if ( temp2<min_dis2 )
                            min_dis2=temp2;
                            min_dis_cluster2=c2;
                        end
                    end
                    if (min_dis2>do)
                        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis2 * min_dis2 * min_dis2 * min_dis2));
                    else
                        S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis2 * min_dis2));
                    end
                    if(min_dis2>0)
                        S(C2(min_dis_cluster2).id).E = S(C2(min_dis_cluster2).id).E- ( (ERX + EDA)*4000 ); 
%                         PACKETS_TO_CH1(r+1)=n-dead-cluster1+1;
                    end
                    S(i).min_dis2=min_dis2;
                    S(i).min_dis_cluster2=min_dis_cluster2;
                end
                
            elseif ismember(i,S3)
                if(cluster3-1>=1)
                    min_dis3=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                    min_dis_cluster3=1;
                    for c3=1:1:cluster3-1
                        temp3=min(min_dis3,sqrt( (S(i).xd-C3(c3).xd)^2 + (S(i).yd-C3(c3).yd)^2 ) );
                        if ( temp3<min_dis3 )
                            min_dis3=temp3;
                            min_dis_cluster3=c3;
                        end
                    end
                    if (min_dis3>do)
                        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis3 * min_dis3 * min_dis3 * min_dis3));
                    else
                        S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis3 * min_dis3));
                    end
                    if(min_dis3>0)
                        S(C3(min_dis_cluster3).id).E = S(C3(min_dis_cluster3).id).E- ( (ERX + EDA)*4000 ); 
%                         PACKETS_TO_CH1(r+1)=n-dead-cluster1+1;
                    end
                    S(i).min_dis3=min_dis3;
                    S(i).min_dis_cluster3=min_dis_cluster3;
                end
                
            else
                if(cluster4-1>=1)
                    min_dis4=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                    min_dis_cluster4=1;
                    for c4=1:1:cluster4-1
                        temp4=min(min_dis4,sqrt( (S(i).xd-C4(c4).xd)^2 + (S(i).yd-C4(c4).yd)^2 ) );
                        if ( temp4<min_dis4 )
                            min_dis4=temp4;
                            min_dis_cluster4=c4;
                        end
                    end
                    if (min_dis4>do)
                        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis4 * min_dis4 * min_dis4 * min_dis4));
                    else
                        S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis4 * min_dis4));
                    end
                    if(min_dis4>0)
                        S(C4(min_dis_cluster4).id).E = S(C4(min_dis_cluster4).id).E- ( (ERX + EDA)*4000 ); 
%                         PACKETS_TO_CH1(r+1)=n-dead-cluster1+1;
                    end
                    S(i).min_dis4=min_dis4;
                    S(i).min_dis_cluster4=min_dis_cluster4;
                end
                
            end
            
        end
    end
        
 
    
end


% LEACH

%counter for CHs
countCHs5=0;
%counter for CHs per round
rcountCHs5=0;
cluster5=1;

countCHs5;
rcountCHs5=rcountCHs5+countCHs5;
flag_first_dead5=0;

for g=0:1:rmax
    g

  %Operation for epoch
  if(mod(g, round(1/p) )==0)
    for i=1:1:n
        S5(i).G=0;        
    end
  end



%Number of dead nodes
dead5=0;


%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS5=0;
packets_TO_CH5=0;
%counter for bit transmitted to Bases Station and to Cluster Heads 
%per round
PACKETS_TO_CH5(g+1)=0;
PACKETS_TO_BS5(g+1)=0;



for i=1:1:n
    %checking if there is a dead node
    if (S5(i).E<=0)        
        dead5=dead5+1;          
    end
    if S5(i).E>0
        S5(i).type='N';        
    end
end



STATISTICS(g+1).DEAD5=dead5;
DEAD5(g+1)=dead5;


%When the first node dies
if (dead5==1)
    if(flag_first_dead5==0)
        first_dead5=g;
        flag_first_dead5=1;
    end
end

countCHs5=0;
cluster5=1;
for i=1:1:n
   if(S5(i).E>0)
   temp_rand=rand;     
   if ( (S5(i).G)<=0)

 %Election of Cluster Heads
 if(temp_rand<= (p/(1-p*mod(g,round(1/p)))))
            countCHs5=countCHs5+1;
            packets_TO_BS5=packets_TO_BS5+1;
            PACKETS_TO_BS5(g+1)=packets_TO_BS5;
            
            S5(i).type='C';
            S5(i).G=round(1/p)-1;
            C5(cluster5).xd=S5(i).xd;
            C5(cluster5).yd=S5(i).yd;           
            
            distance5=sqrt( (S5(i).xd-(S5(n+1).xd) )^2 + (S5(i).yd-(S5(n+1).yd) )^2 );
            C5(cluster5).distance=distance5;
            C5(cluster5).id=i;
            X5(cluster5)=S5(i).xd;
            Y5(cluster5)=S5(i).yd;
            cluster5=cluster5+1;
            
            %Calculation of Energy dissipated
            distance5;
            if (distance5>do)
                S5(i).E=S5(i).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance5*distance5*distance5*distance5 )); 
            end
            if (distance5<=do)
                S5(i).E=S5(i).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance5 * distance5 )); 
            end
        end     
    
    end
  end 
end

STATISTICS(g+1).CLUSTERHEADS5=cluster5-1;
CLUSTERHS5(g+1)=cluster5-1;

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S5(i).type=='N' && S5(i).E>0 )
     if(cluster5-1>=1)
       min_dis5=sqrt( (S5(i).xd-S5(n+1).xd)^2 + (S5(i).yd-S5(n+1).yd)^2 );
       min_dis_cluster5=1;
       for c5=1:1:cluster5-1
           temp5=min(min_dis5,sqrt( (S5(i).xd-C5(c5).xd)^2 + (S5(i).yd-C5(c5).yd)^2 ) );
           if ( temp5<min_dis5 )
               min_dis5=temp5;
               min_dis_cluster5=c5;
           end
       end
       
       %Energy dissipated by associated Cluster Head
            min_dis5;
            if (min_dis5>do)
                S5(i).E=S5(i).E- ( ETX*(4000) + Emp*4000*( min_dis5 * min_dis5 * min_dis5 * min_dis5)); 
            end
            if (min_dis5<=do)
                S5(i).E=S5(i).E- ( ETX*(4000) + Efs*4000*( min_dis5 * min_dis5)); 
            end
        %Energy dissipated
        if(min_dis5>0)
          S5(C5(min_dis_cluster5).id).E = S5(C5(min_dis_cluster5).id).E- ( (ERX + EDA)*4000 ); 
         PACKETS_TO_CH5(g+1)=n-dead5-cluster5+1; 
        end

       S5(i).min_dis=min_dis5;
       S5(i).min_dis_cluster=min_dis_cluster5;
           
   end
 end
end


countCHs5;
rcountCHs5=rcountCHs5+countCHs5;


end


























figure(11);
x1=1:1:r;
y1=1:1:r;
x2=1:1:g;
y2=1:1:g;

for i=1:rmax;
    x1(i)=i;
    y1(i) = n-STATISTICS(i).DEAD;  
    x2(i)=i;
    y2(i) = n - STATISTICS(i).DEAD5;
end

plot(x1,y1,'r--',x2,y2,'g--');

