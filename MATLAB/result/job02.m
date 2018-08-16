clear;
for time=1:1:2
   

time
xm=500;
ym=500;
sink.x=xm*0.5;
sink.y=ym*0.5;
n=1000;
p1=0.05;
p2=0.01;
packetLength = 4000;%数据包长度
ctrPacketLength = 100;%控制包长度
Eo=0.5;
ETX=50*0.000000001;
ERX=50*0.000000001;
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
EDA=5*0.000000001;
rmax=3000;
do=sqrt(Efs/Emp);

C=cell(1,30);%存fv
B=cell(1,30);
index=cell(1,30);
Rank=cell(1,30);
k=2;
b1=[];
b2=[];
cluster=2;

%选取节点
%figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;%坐标
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;    
    S(i).id=i;%编号
    S(i).G=0;
    S(i).E=Eo;
    S(i).ENERGY=0;
    S(i).type='N';
    %subplot(2,2,1);
    %plot(S(i).xd,S(i).yd,'ko');
    %hold on; 
    %test1=num2str(i);
    %test1=[' ',test1];
    %text(XR(i),YR(i),test1)
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
%subplot(2,2,1);
%plot(S(n+1).xd,S(n+1).yd,'rx');
%hold off;

S1=S;
S2=S;
for i=1:n
    S2(i).Dead=0;
end

%执行单一leach
countCHs1=0;%counter for CHs
rcountCHs1=0;%counter for CHs per round
cluster1=1;

countCHs1;
rcountCHs1=rcountCHs1+countCHs1;
flag_first_dead1=0;

for r1=0:1:rmax %主循环,每次1轮
    time
    r1
    %Operation for epoch
    if(mod(r1, round(1/p1) )==0)%round(1/p)取整，四舍五入
        for i=1:1:n
            S1(i).G=0;
            S1(i).cl=0;
        end
    end
    
    dead1=0;%Number of dead nodes
    
    %counter for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS1=0;
    packets_TO_CH1=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads per round
    PACKETS_TO_CH1(r1+1)=0;
    PACKETS_TO_BS1(r1+1)=0;
    
    for i=1:1:n
        %checking if there is a dead node
        if (S1(i).E<=0)    
            dead1=dead1+1;
        end
        if S1(i).E>0
            S1(i).type='N';   
        end
    end
    
    STATISTICS1(time,r1+1).DEAD=dead1;
    DEAD1(r1+1)=dead1;
    
    %When the first node dies
    if (dead1==1)
        if(flag_first_dead1==0)
            first_dead1=r1;
            flag_first_dead1=1;
        end
    end
    
    %计算消耗能量
    Remain_E1=0;
    for i=1:n
        Remain_E1=Remain_E1+S1(i).E;
    end
    
    STATISTICS.Remain_E1(time,r1+1)=Remain_E1;   

    
    countCHs1=0;
    cluster1=1;
for i=1:1:n
   if(S1(i).E>0)
     temp_rand1=rand;    
     if ( (S1(i).G)<=0) %如果该节点在候选集合中
        %Election of Cluster Heads
        if( temp_rand1 <= (p1/(1-p1*mod(r1,round(1/p1)))))
            countCHs1 = countCHs1+1;           
            S1(i).type = 'C';
            S1(i).G = round(1/p1)-1;
            CL1(cluster1).xd = S1(i).xd;
            CL1(cluster1).yd = S1(i).yd;            
           
            distance1=sqrt( (S1(i).xd-(S(n+1).xd) )^2 + (S1(i).yd-(S(n+1).yd) )^2 );%到sink的距离
            CL1(cluster1).distance = distance1;
            CL1(cluster1).id = i;
            X1(cluster1)=S1(i).xd;
            Y1(cluster1)=S1(i).yd;
            cluster1=cluster1+1;
            %%广播自成为簇头
            distanceBroad1 = sqrt(xm*xm+ym*ym);
            if (distanceBroad1 > do)
                S1(i).E = S1(i).E- ( ETX * ctrPacketLength + Emp* ctrPacketLength*( distanceBroad1*distanceBroad1*distanceBroad1*distanceBroad1 ));%广播自成为簇头
            else
                S1(i).E = S1(i).E- ( ETX * ctrPacketLength + Efs * ctrPacketLength*( distanceBroad1*distanceBroad1));
            end
            %Calculation of Energy dissipated 簇头自己发送数据包能量消耗
            distance1;
            if (distance1 > do)
                 S1(i).E = S1(i).E- ( (ETX+EDA)*(packetLength) + Emp * packetLength*( distance1*distance1*distance1*distance1 ));
            else
                 S1(i).E = S1(i).E- ( (ETX+EDA)*(packetLength) + Efs * packetLength*( distance1 * distance1 ));
            end
            packets_TO_BS1 = packets_TO_BS1+1;
            PACKETS_TO_BS1(r1+1) = packets_TO_BS1;
        end    
   
     end
   end
end

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S1(i).type=='N' && S1(i).E>0 ) %普通节点   
     min_dis1 = inf;
     if(cluster1 -1 >= 1)%如果有簇头存在
         min_dis_cluster1 = 1;
         %加入最近的簇头
         for c = 1:1:cluster1 - 1 %簇头数量一共是cluster - 1            
            temp1 = sqrt( (S1(i).xd - CL1(c).xd)^2 + (S1(i).yd - CL1(c).yd)^2 );
            if ( temp1 < min_dis1 )
                min_dis1 = temp1;
                min_dis_cluster1 = c;
            end
            %接收簇头发来的广播的消耗
            S1(i).E = S1(i).E - ETX * ctrPacketLength;
         end
      
         %Energy dissipated by associated Cluster Head普通节点发送数据包到簇头消耗,和加入消息
         min_dis1;
         if (min_dis1 > do)
             S1(i).E = S1(i).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1)); %向簇头发送加入控制消息
             S1(i).E = S1(i).E - ( ETX*(packetLength) + Emp*packetLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1)); %向簇头数据包
         else
            S1(i).E = S1(i).E - ( ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis1 * min_dis1)); %向簇头发送加入控制消息
            S1(i).E = S1(i).E - ( ETX*(packetLength) + Efs*packetLength*( min_dis1 * min_dis1)); %向簇头数据包
         end
         S1(i).E = S1(i).E - ETX*(ctrPacketLength); %接收簇头确认加入控制消息
            
         %Energy dissipated %簇头接收簇成员数据包消耗能量,接收加入消息和和确认加入消息
         if(min_dis1 > 0)
            S1(CL1(min_dis_cluster1).id).E = S1(CL1(min_dis_cluster1).id).E - ( (ERX + EDA)*packetLength ); %接受簇成员发来的数据包
            S1(CL1(min_dis_cluster1).id).E = S1(CL1(min_dis_cluster1).id).E - ERX *ctrPacketLength ; %接收加入消息
            if (min_dis1 > do)%簇头向簇成员发送确认加入的消息
                S1(CL1(min_dis_cluster1).id).E = S1(CL1(min_dis_cluster1).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1));
            else
                S1(CL1(min_dis_cluster1).id).E = S1(CL1(min_dis_cluster1).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis1 * min_dis1));
            end
           PACKETS_TO_CH1(r1+1) = n - dead1 - cluster1 + 1; %所有的非死亡的普通节点都发送数据包
         end
      
         S1(i).min_dis = min_dis1;
         S1(i).min_dis_cluster = min_dis_cluster1;
     else
         min_dis1=sqrt( (S1(i).xd - S(n+1).xd)^2 + (S1(i).yd - S(n+1).yd)^2 );
         if (min_dis1>do)
             S1(i).E=S1(i).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1));
         end
         if (min_dis1<=do)
             S1(i).E=S1(i).E- ( ETX*(packetLength) + Efs*packetLength*( min_dis1 * min_dis1)); 
         end
         packets_TO_BS1=packets_TO_BS1+1;
     end
   end
end

STATISTICS.PACKETS_TO_CH1(time,r1+1)=packets_TO_CH1;
STATISTICS.PACKETS_TO_BS1(time,r1+1)=packets_TO_BS1;

countCHs1;
rcountCHs1 = rcountCHs1 + countCHs1;

end



%分区
for i=1:n
    for j=1:n
        adj(i,j)=sqrt((XR(i)-XR(j))^2 + (YR(i)-YR(j))^2 );
        if adj(i,j)> 50
            adj(i,j)=0;   
        elseif i==j 
            adj(i,j)=0;
        else adj(i,j)=1;
        end
    end
end

fv=fiedlerVector(adj);
A1=find(fv>=0);
A2=find(fv<0);
for i=1:length(adj)
   if ismember(i,A1)
       b1=[b1,i];
   end
end
for i=1:length(adj)
   if ismember(i,A2)
       b2=[b2,i];
   end
end
B{1}=b1;
B{2}=b2;
for i=1:length(B{1})
    for j=1:length(B{1})
        C1(i,j)=sqrt((XR(B{1}(1,i))-XR(B{1}(1,j)))^2 + (YR(B{1}(1,i))-YR(B{1}(1,j)))^2 );
        if C1(i,j)>50
            C1(i,j)=0;   
        elseif i==j 
            C1(i,j)=0;
        else C1(i,j)=1;
        end
    end
end
C{1}=C1;
s1=graphSpectrum(C1);
for i=1:length(B{2})
    for j=1:length(B{2})
        C2(i,j)=sqrt((S(B{2}(1,i)).xd-S(B{2}(1,j)).xd)^2 + (S(B{2}(1,i)).yd-S(B{2}(1,j)).yd)^2 );
        if C2(i,j)>50
            C2(i,j)=0;   
        elseif i==j 
            C2(i,j)=0;
        else C2(i,j)=1;
        end
    end
end
C{2}=C2;
s2=graphSpectrum(C2);
P=[s1,s2];

% figure(1);
% for i=1:length(B{1})
%     subplot(2,2,2);
%     plot(XR(B{1}(1,i)),YR(B{1}(1,i)),'ko');
%     hold on;
% end
% for i=1:length(B{2})
%     subplot(2,2,2);
%     plot(XR(B{2}(1,i)),YR(B{2}(1,i)),'ro');
%     hold on;
% end
% subplot(2,2,2);
% plot(S(n+1).xd,S(n+1).yd,'rx');


while (cluster<4)    
    cluster=cluster+1;
    m=length(P);
    P_min=min(P);
    for e=1:m  
        if P(e)==P_min
            P(e)=inf;
            b1=[];
            b2=[];
            A1=[];
            A2=[];
            fv=fiedlerVector(C{e});
            A1=find(fv>=0);
            A2=find(fv<0);
            for i=1:length(C{e})
                if ismember(i,A1)
                    b1=[b1,B{e}(1,i)];
                end
            end
            for i=1:length(C{e})
                if ismember(i,A2)
                    b2=[b2,B{e}(1,i)];
                end
            end
            C{e}=[];
            B{e}=[];
            k=k+1;
            B{k}=b1;
            C1=[];
            for i=1:length(b1)
                for j=1:length(b1)                    
                    C1(i,j)=sqrt((S(b1(i)).xd-S(b1(j)).xd)^2 + (S(b1(i)).yd-S(b1(j)).yd)^2 );
                    if C1(i,j)>50
                        C1(i,j)=0;
                    elseif i==j
                        C1(i,j)=0;
                    else C1(i,j)=1;
                    end
                end
            end
            C{k}=C1;
            s1=graphSpectrum(C1);
            k=k+1;            
            B{k}=b2;
            C2=[];
            for i=1:length(b2)
                for j=1:length(b2)                    
                    C2(i,j)=sqrt((S(b2(i)).xd-S(b2(j)).xd)^2 + (S(b2(i)).yd-S(b2(j)).yd)^2 );
                    if C2(i,j)>50
                        C2(i,j)=0;
                    elseif i==j
                        C2(i,j)=0;
                    else C2(i,j)=1;
                    end
                end
            end
            C{k}=C2;
            s2=graphSpectrum(C2);
            P=[P,s1,s2];            
        end       
    end
%     figure(1);
%     for i=1:length(B)
%         if ~isempty(B{i})
%             switch i
%                 case {1}
%                     for j=1:length(B{1})
%                         subplot(2,2,cluster);
%                         plot(XR(B{1}(1,j)),YR(B{1}(1,j)),'ko');
%                         hold on;
%                     end
%                 case {2}
%                     for j=1:length(B{2})
%                         subplot(2,2,cluster);
%                         plot(XR(B{2}(1,j)),YR(B{2}(1,j)),'r*');
%                         hold on;
%                     end
%                 case {3}
%                     for j=1:length(B{3})
%                         subplot(2,2,cluster);
%                         plot(XR(B{3}(1,j)),YR(B{3}(1,j)),'g*');
%                         hold on;
%                     end    
%                  case {4}
%                     for j=1:length(B{4})
%                         subplot(2,2,cluster);
%                         plot(XR(B{4}(1,j)),YR(B{4}(1,j)),'b*');
%                         hold on;
%                     end 
%                  case {5}
%                     for j=1:length(B{5})
%                         subplot(2,2,cluster);
%                         plot(XR(B{5}(1,j)),YR(B{5}(1,j)),'co');
%                         hold on;
%                     end 
%                  case {6}
%                     for j=1:length(B{6})
%                         subplot(2,2,cluster);
%                         plot(XR(B{6}(1,j)),YR(B{6}(1,j)),'k*');
%                         hold on;
%                     end   
%                     case {7}
%                     for j=1:length(B{7})
%                         subplot(2,2,cluster);
%                         plot(XR(B{7}(1,j)),YR(B{7}(1,j)),'mo');
%                         hold on;
%                     end
%                     case {8}
%                     for j=1:length(B{8})
%                         subplot(2,2,cluster);
%                         plot(XR(B{8}(1,j)),YR(B{8}(1,j)),'ko');
%                         hold on;
%                     end
%                 case {9}
%                     for j=1:length(B{9})
%                         subplot(2,2,cluster);
%                         plot(XR(B{9}(1,j)),YR(B{9}(1,j)),'ro');
%                         hold on;
%                     end
%                 case {10}
%                     for j=1:length(B{10})
%                         subplot(2,2,cluster);
%                         plot(XR(B{10}(1,j)),YR(B{10}(1,j)),'go');
%                         hold on;
%                     end    
%                  case {11}
%                     for j=1:length(B{11})
%                         subplot(2,2,cluster);
%                         plot(XR(B{11}(1,j)),YR(B{11}(1,j)),'bo');
%                         hold on;
%                     end 
%                  case {12}
%                     for j=1:length(B{12})
%                         subplot(2,2,cluster);
%                         plot(XR(B{12}(1,j)),YR(B{12}(1,j)),'co');
%                         hold on;
%                     end 
%                  case {13}
%                     for j=1:length(B{13})
%                         subplot(2,2,cluster);
%                         plot(XR(B{13}(1,j)),YR(B{13}(1,j)),'yo');
%                         hold on;
%                     end   
%                     case {14}
%                     for j=1:length(B{14})
%                         subplot(2,2,cluster);
%                         plot(XR(B{14}(1,j)),YR(B{14}(1,j)),'mo');
%                         hold on;
%                     end
%             end
%         end
%     end
%     subplot(2,2,cluster);
%     plot(S(n+1).xd,S(n+1).yd,'rx');
end

%选择leader closeness
for i=1:length(C)
    if ~isempty(C{i})        
        CC=zeros(length(C{i}),1);
        for j=1:length(C{i})
            CC(j)=1/sum( simple_dijkstra(C{i},j) ); 
        end
        [Rank1,index1]= sort(CC);
        index{i}=index1;
        Rank{i}=Rank1;
    end
end
for i=1:length(index)
    if ~isempty(index{i})
        %figure(cluster);        
        %plot(S2(B{i}(index{i}(1,1))).xd,S2(B{i}(index{i}(1,1))).yd,'k*');
        S2(B{i}(index{i}(1,1))).E = 1*Eo;
        S2(B{i}(index{i}(1,1))).type='L';
        
    end
end


%分区后执行leach
%counter for CHs
countCHs2=0;
%counter for CHs per round
rcountCHs2=0;
cluster2=1;

countCHs2;
rcountCHs2=rcountCHs2+countCHs2;
flag_first_dead2=0;
dead2=0;

for r2=0:1:rmax
    time
    r2
    if(mod(r2, round(1/p2) )==0)
        for i=1:1:n
            if S2(i).type ~='L'
                S2(i).G=0;
                S2(i).cl=0;
            end
        end
    end
    
    packets_TO_BS2=0;
    packets_TO_CH2=0;
    
    PACKETS_TO_CH2(r2+1)=0;
    PACKETS_TO_BS2(r2+1)=0;
    
    for i=1:n
        if (S2(i).E<=0 && S2(i).Dead==0)
            dead2=dead2+1;
            S2(i).Dead=1;
        end
        if S2(i).type~='L' && S2(i).E>0
            S2(i).type='N';
        end
    end
    
    STATISTICS2(time,r2+1).DEAD=dead2;
    DEAD2(r2+1)=dead2;
    
    if (dead2==1)
        if(flag_first_dead2==0)
            first_dead2=r2;
            flag_first_dead2=1;
        end
    end
    
    %计算消耗能量
    Remain_E2=0;
    for i=1:n
        Remain_E2=Remain_E2+S2(i).E;
    end
   
    STATISTICS.Remain_E2(time,r2+1)=Remain_E2;
    
    
    countCHs2=0;
    cluster2=1;
    
    for i=1:length(B)
        if ~isempty(B{i})
            %计算leader到BS的距离
            distance_Leader_TO_BS = sqrt((S2(B{i}(index{i}(1,1))).xd-S(n+1).xd)^2 + (S2(B{i}(index{i}(1,1))).yd-S(n+1).yd)^2);
            
            for j=1:1:length(B{i})
                if (S2(B{i}(1,j)).E>0 && S2(B{i}(1,j)).type=='N')
                    temp_rand2=rand;
                    if ( (S2(B{i}(1,j)).G)<=0) 
                        if( temp_rand2 <= (p2/(1-p2*mod(r2,round(1/p2)))))
                            countCHs2 = countCHs2+1;
                            S2(B{i}(1,j)).type = 'C';
                            S2(B{i}(1,j)).G = round(1/p2)-1;
                            CL2(cluster2).xd = S2(B{i}(1,j)).xd;
                            CL2(cluster2).yd = S2(B{i}(1,j)).yd;
                            distance2=sqrt( (S2(B{i}(1,j)).xd-(S2(B{i}(index{i}(1,1))).xd) )^2 + (S2(B{i}(1,j)).yd-(S2(B{i}(index{i}(1,1))).yd) )^2 );%到leader的距离
                            CL2(cluster2).distance = distance2;
                            CL2(cluster2).id = B{i}(1,j);
                            X2(cluster2)=S(B{i}(1,j)).xd;
                            Y2(cluster2)=S(B{i}(1,j)).yd;
                            cluster2=cluster2+1;
                            
                            distanceBroad2 = sqrt(xm*xm/4+ym*ym/4); %CH广播
                            if (distanceBroad2 > do)
                                S2(B{i}(1,j)).E = S2(B{i}(1,j)).E- ( ETX * ctrPacketLength + Emp* ctrPacketLength*( distanceBroad2*distanceBroad2*distanceBroad2*distanceBroad2 ));
                            else
                                S2(B{i}(1,j)).E = S2(B{i}(1,j)).E- ( ETX * ctrPacketLength + Efs * ctrPacketLength*( distanceBroad2*distanceBroad2));
                            end
                            
                            distance2;%CH到Leader的距离
                            if (distance2 > do)
                                S2(B{i}(1,j)).E = S2(B{i}(1,j)).E- ( (ETX+EDA)*(packetLength) + Emp * packetLength*( distance2*distance2*distance2*distance2 ));
                            else
                                S2(B{i}(1,j)).E = S2(B{i}(1,j)).E- ( (ETX+EDA)*(packetLength) + Efs * packetLength*( distance2 * distance2 ));
                            end
                            packets_TO_BS2 = packets_TO_BS2+1;
                            PACKETS_TO_BS2(r2+1) = packets_TO_BS2;
                        end
                    end
                end
            end
            
            for j=1:1:length(B{i})
                if ( S2(B{i}(1,j)).type=='N' && S2(B{i}(1,j)).E>0 ) 
                    min_dis2 = inf;
                    if(cluster2 -1 >= 1)%如果有簇头存在
                        min_dis_cluster2 = 1;
                        %加入最近的簇头
                        for c = 1:1:cluster2 - 1 %簇头数量一共是cluster - 1
                            %temp = min(min_dis,sqrt( (S(i).xd - C(c).xd)^2 + (S(i).yd - C(c).yd)^2 ) );
                            temp2 = sqrt( (S2(B{i}(1,j)).xd - CL2(c).xd)^2 + (S2(B{i}(1,j)).yd - CL2(c).yd)^2 );
                            if ( temp2 < min_dis2 )
                                min_dis2 = temp2;
                                min_dis_cluster2 = c;
                            end
                            %接收簇头发来的广播的消耗
                            S2(B{i}(1,j)).E = S2(B{i}(1,j)).E - ETX * ctrPacketLength;
                        end
                        
                        %Energy dissipated by associated Cluster Head普通节点发送数据包到簇头消耗,和加入消息
                        min_dis2;
                        if (min_dis2 > do)
                            S2(B{i}(1,j)).E = S2(B{i}(1,j)).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2)); %向簇头发送加入控制消息
                            S2(B{i}(1,j)).E = S2(B{i}(1,j)).E - ( ETX*(packetLength) + Emp*packetLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2)); %向簇头数据包
                        else
                            S2(B{i}(1,j)).E = S2(B{i}(1,j)).E - ( ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis2 * min_dis2)); %向簇头发送加入控制消息
                            S2(B{i}(1,j)).E = S2(B{i}(1,j)).E - ( ETX*(packetLength) + Efs*packetLength*( min_dis2 * min_dis2)); %向簇头数据包
                        end
                        S2(B{i}(1,j)).E = S2(B{i}(1,j)).E - ETX*(ctrPacketLength); %接收簇头确认加入控制消息
                        
                        %Energy dissipated %簇头接收簇成员数据包消耗能量,接收加入消息和和确认加入消息
                        if(min_dis2 > 0)
                            S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( (ERX + EDA)*packetLength ); %接受簇成员发来的数据包
                            S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ERX *ctrPacketLength ; %接收加入消息
                            if (min_dis2 > do)%簇头向簇成员发送确认加入的消息
                                S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2));
                            else
                                S(CL2(min_dis_cluster2).id).E = S(CL2(min_dis_cluster2).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis2 * min_dis2));
                            end
                            PACKETS_TO_CH2(r2+1) = n - dead2 - cluster2 + 1; %所有的非死亡的普通节点都发送数据包
                            if (distance_Leader_TO_BS>do)
                                S2(B{i}(index{i}(1,1))).E = S2(B{i}(index{i}(1,1))).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( distance_Leader_TO_BS * distance_Leader_TO_BS * distance_Leader_TO_BS * distance_Leader_TO_BS));
                            else
                                S2(B{i}(index{i}(1,1))).E = S2(B{i}(index{i}(1,1))).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( distance_Leader_TO_BS * distance_Leader_TO_BS ));
                            end
                        end
                        S2(B{i}(1,j)).min_dis = min_dis2;
                        S2(B{i}(1,j)).min_dis_cluster = min_dis_cluster2;
                    else
                        min_dis2=sqrt( (S2(B{i}(1,j)).xd - S2(B{i}(index{i}(1,1))).xd)^2 + (S2(B{i}(1,j)).yd - S2(B{i}(index{i}(1,1))).yd)^2 );
                        if (min_dis2>do)
                            S2(B{i}(1,j)).E=S2(B{i}(1,j)).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis2 * min_dis2 * min_dis2 * min_dis2));
                        end
                        if (min_dis2<=do)
                            S2(B{i}(1,j)).E=S2(B{i}(1,j)).E- ( ETX*(packetLength) + Efs*packetLength*( min_dis2 * min_dis2));
                        end
                        packets_TO_BS2=packets_TO_BS2+1;
                    end
                end
            end
            
            STATISTICS.PACKETS_TO_CH2(time,r2+1)=packets_TO_CH2;
            STATISTICS.PACKETS_TO_BS2(time,r2+1)=packets_TO_BS2;
            countCHs2;
            rcountCHs2 = rcountCHs2 + countCHs2;     
            
            
        end
    end
end


%LEACH-C
p3=0.05;
A3=zeros(1,rmax+1);
D3=zeros(1,rmax+1);
T3=zeros(1,rmax+1);
MIN3 = zeros(1,rmax+1);
MAX3 = zeros(1,rmax+1);

%产生随机节点
for i=1:n
    S3(i).xd=XR(i);
    S3(i).yd=YR(i);
    S3(i).G=0;
    S3(i).type='N';
    S3(i).E=Eo;
    S3(i).ENERGY=0;
    
end

mind3=sqrt( (S3(1).xd-(S(n+1).xd) )^2 + (S3(1).yd-(S(n+1).yd) )^2 );
maxd3=sqrt( (S3(1).xd-(S(n+1).xd) )^2 + (S3(1).yd-(S(n+1).yd) )^2 );
min_id3=1;
max_id3=1;
for i=2:1:n
    tempd3=sqrt( (S3(i).xd-(S(n+1).xd) )^2 + (S3(i).yd-(S(n+1).yd) )^2 );
    if(tempd3>maxd3)
    maxd3=tempd3;
    max_id3=i;
    end
    if(tempd3<mind3)
    mind3=tempd3;
    min_id3=i;
    end
end

countCHs3=0;
rcountCHs3=0;
cluster3=1;

countCHs3;
rcountCHs3=rcountCHs3+countCHs3;
flag_first_dead3=0;

for r3=0:1:rmax
    time
    r3
    
    E3=0;
    if(mod(r3, round(1/p3) )==0)
        for i=1:1:n
            S3(i).G=0;
            S3(i).cl=0;
        end
    end
    
    dead3=0;
    dead_a3=0;
    dead_n3=0;
    alive3=n;
    
    packets_TO_BS3=0;
    packets_TO_CH3=0;
    PACKETS_TO_CH3(r3+1)=0;
    PACKETS_TO_BS3(r3+1)=0;
    
    if(S3(min_id3).E>0)
        MIN3(r3+1)=S3(min_id3).E; 
    end
    if(S3(max_id3).E>0)
        MAX3(r3+1)=S3(max_id3).E;
    end
    

    for i=1:1:n
        if (S3(i).E<=0)
            S3(i).E=0;
            dead3=dead3+1;
            alive3=alive3-1;
            if(S3(i).ENERGY==1)
                dead_a3=dead_a3+1;
            end
            if(S3(i).ENERGY==0)
                dead_n3=dead_n3+1;
            end
        end
        if S3(i).E>0
            S3(i).type='N';            
        end
    end
    STATISTICS3(time,r3+1).DEAD=dead3;
    DEAD3(r3+1)=dead3;
    D3(r3+1)=dead3;
    DEAD_N3(r3+1)=dead_n3;
    DEAD_A3(r3+1)=dead_a3;
    ALIVE3(r3+1)=alive3;
    
    if (dead3==1)
        if(flag_first_dead3==0)
            first_dead3=r3;
            flag_first_dead3=1;
        end
    end
    
    %计算消耗能量
    Remain_E3=0;
    for i=1:n
        Remain_E3=Remain_E3+S3(i).E;
    end
    
    STATISTICS.Remain_E3(time,r3+1)=Remain_E3;
   
    
    
    countCHs3=0;
    cluster3=1;
    
    for i=1:n;
        if(S3(i).E>0)
            E3=E3+S3(i).E;
        end
    end
    
    for i=1:1:n
        if(S3(i).E>0)
            temp_rand3=rand; 
            
            if ( (S3(i).G)<=0)
                if(temp_rand3<= ((S3(i).E*alive3*p3)/E3)) 
                    countCHs3=countCHs3+1;
                    packets_TO_BS3=packets_TO_BS3+1;
                    PACKETS_TO_BS3(r3+1)=packets_TO_BS3;
                    S3(i).type='C';
                    S3(i).G=round(1/p3)-1; 
                    CL3(cluster3).xd=S3(i).xd;
                    CL3(cluster3).yd=S3(i).yd;
                    distance3=sqrt( (S3(i).xd-(S(n+1).xd) )^2 + (S3(i).yd-(S(n+1).yd) )^2 );                    
                    CL3(cluster3).distance=distance3;
                    CL3(cluster3).id=i;
                    X3(cluster3)=S3(i).xd;
                    Y3(cluster3)=S3(i).yd;
                    cluster3=cluster3+1;
                    
                    
                    distance3;
                    if (distance3>do) 
                        S3(i).E=S3(i).E- ( (ETX+EDA)*(packetLength) + Emp*packetLength*( distance3*distance3*distance3*distance3 )); 
                    else
                        S3(i).E=S3(i).E- ( (ETX+EDA)*(packetLength) + Emp*packetLength*( distance3*distance3 ));
                    end
                end
            end
        end
    end
    
    %STATISTICS(r3+1).CLUSTERHEADS=cluster3-1;
    %CLUSTERHS(r3+1)=cluster3-1;
    
    for i=1:1:n
        if ( S3(i).type=='N' && S3(i).E>0 ) 
            if(cluster3-1>=1)
                min_dis3=inf;
                min_dis_cluster3=1;
                for c=1:1:cluster3-1
                    temp3=min(min_dis3,sqrt( (S3(i).xd-CL3(c).xd)^2 + (S3(i).yd-CL3(c).yd)^2 ) );
                    if ( temp3<min_dis3 )
                        min_dis3=temp3;
                        min_dis_cluster3=c;
                    end
                end
                min_dis3; 
                if(min_dis3>0)
                    if (min_dis3>do)
                        S3(i).E=S3(i).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis3 * min_dis3 * min_dis3 * min_dis3)); 
                    else
                        S3(i).E=S3(i).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis3 * min_dis3 ));
                    end
                    dd3=CL3(min_dis_cluster3).distance;
                    if(dd3>do)
                        S3(CL3(min_dis_cluster3).id).E = S3(CL3(min_dis_cluster3).id).E - ( (ERX + EDA)*packetLength + Emp*packetLength*( dd3*dd3*dd3*dd3 ) );
                    else
                        S3(CL3(min_dis_cluster3).id).E = S3(CL3(min_dis_cluster3).id).E - ( (ERX + EDA)*packetLength + Emp*packetLength*( dd3*dd3 ) );
                    end
                    PACKETS_TO_CH3(r3+1)=n-dead3-cluster3+1;
                    PACKETS_TO_BS3(r3+1)=n-dead3-cluster3+1;
                end
                packets_TO_BS3=packets_TO_BS3+1;
                S3(i).min_dis=min_dis3;
                S3(i).min_dis_cluster=min_dis_cluster3;
            end
            if(cluster3-1==0)
                ddd0 = sqrt((S3(i).xd-S(n+1).xd)^2 + (S3(i).yd-S(n+1).yd)^2);
                if (ddd0>do)
                    S3(i).E=S3(i).E- ( ETX*(packetLength) + Emp*packetLength*( ddd0 * ddd0 * ddd0 * ddd0));
                else
                    S3(i).E=S3(i).E- ( ETX*(packetLength) + Emp*packetLength*( ddd0 * ddd0 ));
                end
                packets_TO_BS3=packets_TO_BS3+1;
            end
        end
    end
    for i=1:1:n
        if(S3(i).E>=0)
            A3(1,r3+1)=A3(1,r3+1)+S3(i).E;
        end
    end
    T3(1,r3+1)=T3(1,r3+1)+packets_TO_BS3;
    countCHs3;
    rcountCHs3=rcountCHs3+countCHs3;
    STATISTICS.PACKETS_TO_CH3(time,r3+1)=packets_TO_CH3;
    STATISTICS.PACKETS_TO_BS3(time,r3+1)=packets_TO_BS3;
end




%SEP
p5=0.1;
m5=0.1;
a5=1;
for i=1:1:n
    S5(i).xd=XR(i);
    S5(i).yd=YR(i);
    S5(i).G=0;
    S5(i).type='N';
    temp_rnd0=i;
    if (temp_rnd0>=m5*n+1) 
        S5(i).E=Eo;
        S5(i).ENERGY=0;
    end
    if (temp_rnd0<m5*n+1) 
        S5(i).E=Eo*(1+a5);
        S5(i).ENERGY=1;
    end
end
countCHs5=0;
rcountCHs5=0;
cluster5=1;

countCHs5;
rcountCHs5=rcountCHs5+countCHs5;
flag_first_dead5=0;


for r5=0:1:rmax
    time
    r5
    
    pnrm=( p5/ (1+a5*m5) );
    padv= ( p5*(1+a5)/(1+a5*m5) );
    if(mod(r5, round(1/pnrm) )==0)
        for i=1:1:n
            S5(i).G=0;
            S5(i).cl=0;
        end
    end
    
    if(mod(r5, round(1/padv) )==0)
        for i=1:1:n
            if(S5(i).ENERGY==1)
                S5(i).G=0;
                S5(i).cl=0;
            end
        end
    end
    
    
    dead5=0;
    dead_a5=0;
    dead_n5=0;
    
    packets_TO_BS5=0;
    packets_TO_CH5=0;
    
    PACKETS_TO_CH5(r5+1)=0;
    PACKETS_TO_BS5(r5+1)=0;
    
    for i=1:1:n
        if (S5(i).E<=0)
            dead5=dead5+1;
            if(S5(i).ENERGY==1)
                dead_a5=dead_a5+1;
            end
            if(S5(i).ENERGY==0)
                dead_n5=dead_n5+1;
            end
        end
        if S5(i).E>0
            S5(i).type='N';
        end
    end
    
    STATISTICS5(time,r5+1).DEAD=dead5;
    DEAD5(r5+1)=dead5;
    DEAD_N5(r5+1)=dead_n5;
    DEAD_A5(r5+1)=dead_a5;
    
    if (dead5==1)
        if(flag_first_dead5==0)
            first_dead5=r5;
            flag_first_dead5=1;
        end
    end
    
    %计算消耗能量
    Remain_E5=0;
    for i=1:n
        Remain_E5=Remain_E5+S5(i).E;
    end
    
    STATISTICS.Remain_E5(time,r5+1)=Remain_E5;
    
    
    
    
    countCHs5=0;
    cluster5=1;
    for i=1:1:n
        if(S5(i).E>0)
            temp_rand5=rand;
            if ( (S5(i).G)<=0)
                if( ( S5(i).ENERGY==0 && ( temp_rand5 <= ( pnrm / ( 1 - pnrm * mod(r5,round(1/pnrm)) )) ) )  )
                    countCHs5=countCHs5+1;
                    packets_TO_BS5=packets_TO_BS5+1;
                    PACKETS_TO_BS5(r5+1)=packets_TO_BS5;
                    S5(i).type='C';
                    S5(i).G=100;
                    CL5(cluster5).xd=S5(i).xd;
                    CL5(cluster5).yd=S5(i).yd;
                    
                    distance5=sqrt( (S5(i).xd-(S(n+1).xd) )^2 + (S5(i).yd-(S(n+1).yd) )^2 );
                    CL5(cluster5).distance=distance5;
                    CL5(cluster5).id=i;
                    X5(cluster5)=S5(i).xd;
                    Y5(cluster5)=S5(i).yd;
                    cluster5=cluster5+1;
                    
                    distance5;
                    if (distance5>do)
                        S5(i).E=S5(i).E- ( (ETX+EDA)*(packetLength) + Emp*packetLength*( distance5*distance5*distance5*distance5 ));
                    else
                        S5(i).E=S5(i).E- ( (ETX+EDA)*(packetLength) + Emp*packetLength*( distance5*distance5 ));
                    end
                end
                if( ( S5(i).ENERGY==1 && ( temp_rand5 <= ( padv / ( 1 - padv * mod(r5,round(1/padv)) )) ) )  )
                    countCHs5=countCHs5+1;
                    packets_TO_BS5=packets_TO_BS5+1;
                    PACKETS_TO_BS5(r5+1)=packets_TO_BS5;
                    
                    S5(i).type='C';
                    S5(i).G=100;
                    CL5(cluster5).xd=S5(i).xd;
                    CL5(cluster5).yd=S5(i).yd;
                    
                    distance5=sqrt( (S5(i).xd-(S(n+1).xd) )^2 + (S5(i).yd-(S(n+1).yd) )^2 );
                    CL5(cluster5).distance=distance5;
                    CL5(cluster5).id=i;
                    X5(cluster5)=S5(i).xd;
                    Y5(cluster5)=S5(i).yd;
                    cluster5=cluster5+1;
                    
                    distance5;
                    if (distance5>do)
                        S5(i).E=S5(i).E- ( (ETX+EDA)*(packetLength) + Emp*packetLength*( distance5*distance5*distance5*distance5 ));
                    else
                        S5(i).E=S5(i).E- ( (ETX+EDA)*(packetLength) + Emp*packetLength*( distance5*distance5 ));
                    end
                end
            end
        end
    end
    %STATISTICS(r5+1).CLUSTERHEADS=cluster5-1;
    %CLUSTERHS(r5+1)=cluster5-1;
    
    for i=1:1:n
        if ( S5(i).type=='N' && S5(i).E>0 )
            if(cluster5-1>=1)
                min_dis5=sqrt( (S5(i).xd-S(n+1).xd)^2 + (S5(i).yd-S(n+1).yd)^2 );
                min_dis_cluster5=1;
                for c=1:1:cluster5-1
                    temp5=min(min_dis5,sqrt( (S5(i).xd-CL5(c).xd)^2 + (S5(i).yd-CL5(c).yd)^2 ) );
                    if ( temp5<min_dis5 )
                        min_dis5=temp5;
                        min_dis_cluster5=c;
                    end
                end
                 min_dis5;
                 if (min_dis5>do)
                     S5(i).E=S5(i).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis5 * min_dis5 * min_dis5 * min_dis5)); 
                 else
                     S5(i).E=S5(i).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis5 * min_dis5 ));
                 end
                 if(min_dis5>0)
                     S5(CL5(min_dis_cluster5).id).E = S5(CL5(min_dis_cluster5).id).E- ( (ERX + EDA)*packetLength );
                     PACKETS_TO_CH5(r5+1)=n-dead5-cluster5+1; 
                 end
                 S5(i).min_dis=min_dis5;
                 S5(i).min_dis_cluster=min_dis_cluster5;
            end
        end
    end
    countCHs5;
    rcountCHs5=rcountCHs5+countCHs5;
    STATISTICS.PACKETS_TO_CH5(time,r5+1)=packets_TO_CH5;
    STATISTICS.PACKETS_TO_BS5(time,r5+1)=packets_TO_BS5;

    
end

end


figure(11);
x11=1:1:r1;
y11=1:1:r1;
x21=1:1:r2;
y21=1:1:r2;
x31=1:1:r3;
y31=1:1:r3;
x51=1:1:r5;
y51=1:1:r5;
Dead_Number1=mean(cell2mat(struct2cell(STATISTICS1)));
Dead_Number2=mean(cell2mat(struct2cell(STATISTICS2)));
Dead_Number3=mean(cell2mat(struct2cell(STATISTICS3)));
Dead_Number5=mean(cell2mat(struct2cell(STATISTICS5)));
for i=1:rmax
    x11(i)=i;
    y11(i) = n - Dead_Number1(i);
    x21(i)=i;
    y21(i) = n - Dead_Number2(i);
    x31(i)=i;
    y31(i) = n - Dead_Number3(i);    
    x51(i)=i;
    y51(i) = n - Dead_Number5(i);
end
subplot(2,1,1);
plot(x11,y11,'r--',x21,y21,'g--',x31,y31,'b--',x51,y51,'k--');
title('图一：存活节点数统计图');
xlabel('工作轮数（/轮）');
ylabel('存活节点数（/个）');
legend('LEACH','NEW','LEACH-C','SEP');




figure(11);
x12=1:1:r1;
y12=1:1:r1;
x22=1:1:r2;
y22=1:1:r2;
x32=1:1:r3;
y32=1:1:r3;
x52=1:1:r5;
y52=1:1:r5;


for i=1:rmax
    x12(i)=i;
    y12(i) = mean(STATISTICS.Remain_E1(:,i));
    x22(i)=i;
    y22(i) = mean(STATISTICS.Remain_E2(:,i));
    x32(i)=i;
    y32(i) = mean(STATISTICS.Remain_E3(:,i));
    x52(i)=i;
    y52(i) = mean(STATISTICS.Remain_E5(:,i));
end
subplot(2,1,2);
plot(x12,y12,'r--',x22,y22,'g--',x32,y32,'b--',x52,y52,'k--');
title('图一：能量消耗统计图');
xlabel('工作轮数（/轮）');
ylabel('剩余能量（/J）');
legend('LEACH','NEW','LEACH-C','SEP');


figure(12);
x13=1:1:r1;
y13=1:1:r1;
x23=1:1:r2;
y23=1:1:r2;
x33=1:1:r3;
y33=1:1:r3;
x53=1:1:r5;
y53=1:1:r5;

for i=1:rmax
    x13(i)=i;
    y13(i) = mean(STATISTICS.PACKETS_TO_BS1(:,i));
    x23(i)=i;
    y23(i) = mean(STATISTICS.PACKETS_TO_BS2(:,i));
    x33(i)=i;
    y33(i) = mean(STATISTICS.PACKETS_TO_BS3(:,i));    
    x53(i)=i;
    y53(i) = mean(STATISTICS.PACKETS_TO_BS5(:,i));
end
subplot(4,1,1);
plot(x13,y13,'r--');
title('A：发送数据统计图');
xlabel('工作轮数（/轮）');
ylabel('数据包数（/个）');
legend('LEACH');

subplot(4,1,2);
plot(x23,y23,'g--');
title('B：发送数据统计图');
xlabel('工作轮数（/轮）');
ylabel('数据包数（/个）');
legend('NEW');

subplot(4,1,3);
plot(x33,y33,'b--');
title('C：发送数据统计图');
xlabel('工作轮数（/轮）');
ylabel('数据包数（/个）');
legend('LEACH-C');

subplot(4,1,4);
plot(x53,y53,'k--');
title('D：发送数据统计图');
xlabel('工作轮数（/轮）');
ylabel('数据包数（/个）');
legend('SEP');