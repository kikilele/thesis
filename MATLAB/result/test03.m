clear;
xm=100;
ym=100;
sink.x=xm*0.5;
sink.y=ym*0.5;
n=200;
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
rmax=6000;
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
figure(1);
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
    plot(S(i).xd,S(i).yd,'ko');
    hold on; 
    %test1=num2str(i);
    %test1=[' ',test1];
    %text(XR(i),YR(i),test1)
end
S(n+1).xd=sink.x;
S(n+1).yd=sink.y;
plot(S(n+1).xd,S(n+1).yd,'rx');
hold off;

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
for r=0:1:rmax %主循环,每次1轮
    r;
    %Operation for epoch
    if(mod(r, round(1/p1) )==0)%round(1/p)取整，四舍五入
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
        end
    end
    
    dead1=0;%Number of dead nodes
    
    %counter for bit transmitted to Bases Station and to Cluster Heads
    packets_TO_BS1=0;
    packets_TO_CH1=0;
    %counter for bit transmitted to Bases Station and to Cluster Heads per round
    PACKETS_TO_CH1(r+1)=0;
    PACKETS_TO_BS1(r+1)=0;
    
    for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)    
            dead1=dead1+1;
        end
        if S(i).E>0
            S(i).type='N';   
        end
    end
    
    STATISTICS1(r+1).DEAD=dead1;
    DEAD1(r+1)=dead1;
    
    %When the first node dies
    if (dead1==1)
        if(flag_first_dead1==0)
            first_dead1=r;
            flag_first_dead1=1;
        end
    end
    
    countCHs1=0;
    cluster1=1;
for i=1:1:n
   if(S(i).E>0)
     temp_rand1=rand;    
     if ( (S(i).G)<=0) %如果该节点在候选集合中
        %Election of Cluster Heads
        if( temp_rand1 <= (p1/(1-p1*mod(r,round(1/p1)))))
            countCHs1 = countCHs1+1;           
            S(i).type = 'C';
            S(i).G = round(1/p1)-1;
            CL1(cluster1).xd = S(i).xd;
            CL1(cluster1).yd = S(i).yd;
            %plot(S(i).xd,S(i).yd,'g*');
           
            distance1=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );%到sink的距离
            CL1(cluster1).distance = distance1;
            CL1(cluster1).id = i;
            X1(cluster1)=S(i).xd;
            Y1(cluster1)=S(i).yd;
            cluster1=cluster1+1;
            %%广播自成为簇头
            distanceBroad1 = sqrt(xm*xm+ym*ym);
            if (distanceBroad1 > do)
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Emp* ctrPacketLength*( distanceBroad1*distanceBroad1*distanceBroad1*distanceBroad1 ));%广播自成为簇头
            else
                S(i).E = S(i).E- ( ETX * ctrPacketLength + Efs * ctrPacketLength*( distanceBroad1*distanceBroad1));
            end
            %Calculation of Energy dissipated 簇头自己发送数据包能量消耗
            distance1;
            if (distance1 > do)
                 S(i).E = S(i).E- ( (ETX+EDA)*(packetLength) + Emp * packetLength*( distance1*distance1*distance1*distance1 ));
            else
                 S(i).E = S(i).E- ( (ETX+EDA)*(packetLength) + Efs * packetLength*( distance1 * distance1 ));
            end
            packets_TO_BS1 = packets_TO_BS1+1;
            PACKETS_TO_BS1(r+1) = packets_TO_BS1;
        end    
   
     end
   end
end

%Election of Associated Cluster Head for Normal Nodes
for i=1:1:n
   if ( S(i).type=='N' && S(i).E>0 ) %普通节点   
     min_dis1 = inf;
     if(cluster1 -1 >= 1)%如果有簇头存在
         min_dis_cluster1 = 1;
         %加入最近的簇头
         for c = 1:1:cluster1 - 1 %簇头数量一共是cluster - 1            
            temp1 = sqrt( (S(i).xd - CL1(c).xd)^2 + (S(i).yd - CL1(c).yd)^2 );
            if ( temp1 < min_dis1 )
                min_dis1 = temp1;
                min_dis_cluster1 = c;
            end
            %接收簇头发来的广播的消耗
            S(i).E = S(i).E - ETX * ctrPacketLength;
         end
      
         %Energy dissipated by associated Cluster Head普通节点发送数据包到簇头消耗,和加入消息
         min_dis1;
         if (min_dis1 > do)
             S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1)); %向簇头发送加入控制消息
             S(i).E = S(i).E - ( ETX*(packetLength) + Emp*packetLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1)); %向簇头数据包
         else
            S(i).E = S(i).E - ( ETX*(ctrPacketLength) + Efs*ctrPacketLength*( min_dis1 * min_dis1)); %向簇头发送加入控制消息
            S(i).E = S(i).E - ( ETX*(packetLength) + Efs*packetLength*( min_dis1 * min_dis1)); %向簇头数据包
         end
         S(i).E = S(i).E - ETX*(ctrPacketLength); %接收簇头确认加入控制消息
            
         %Energy dissipated %簇头接收簇成员数据包消耗能量,接收加入消息和和确认加入消息
         if(min_dis1 > 0)
            S(CL1(min_dis_cluster1).id).E = S(CL1(min_dis_cluster1).id).E - ( (ERX + EDA)*packetLength ); %接受簇成员发来的数据包
            S(CL1(min_dis_cluster1).id).E = S(CL1(min_dis_cluster1).id).E - ERX *ctrPacketLength ; %接收加入消息
            if (min_dis1 > do)%簇头向簇成员发送确认加入的消息
                S(CL1(min_dis_cluster1).id).E = S(CL1(min_dis_cluster1).id).E - ( ETX*(ctrPacketLength) + Emp * ctrPacketLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1));
            else
                S(CL1(min_dis_cluster1).id).E = S(CL1(min_dis_cluster1).id).E - ( ETX*(ctrPacketLength) + Efs * ctrPacketLength*( min_dis1 * min_dis1));
            end
           PACKETS_TO_CH1(r+1) = n - dead1 - cluster1 + 1; %所有的非死亡的普通节点都发送数据包
         end
      
         S(i).min_dis = min_dis1;
         S(i).min_dis_cluster = min_dis_cluster1;
     else
         min_dis1=sqrt( (S(i).xd - S(n+1).xd)^2 + (S(i).yd - S(n+1).yd)^2 );
         if (min_dis1>do)
             S(i).E=S(i).E- ( ETX*(packetLength) + Emp*packetLength*( min_dis1 * min_dis1 * min_dis1 * min_dis1));
         end
         if (min_dis1<=do)
             S(i).E=S(i).E- ( ETX*(packetLength) + Efs*packetLength*( min_dis1 * min_dis1)); 
         end
         packets_TO_BS1=packets_TO_BS1+1;
     end
   end
end

STATISTICS.PACKETS_TO_CH1(r+1)=packets_TO_CH1;
STATISTICS.PACKETS_TO_BS1(r+1)=packets_TO_BS1;

countCHs1;
rcountCHs1 = rcountCHs1 + countCHs1;

end



%分区
for i=1:n
    for j=1:n
        adj(i,j)=sqrt((XR(i)-XR(j))^2 + (YR(i)-YR(j))^2 );
        if adj(i,j)> 30
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
        if C1(i,j)>30
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
        if C2(i,j)>30
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

figure(2);
for i=1:length(B{1})
    plot(XR(B{1}(1,i)),YR(B{1}(1,i)),'ko');
    hold on;
end
for i=1:length(B{2})
    plot(XR(B{2}(1,i)),YR(B{2}(1,i)),'ro');
    hold on;
end
plot(S(n+1).xd,S(n+1).yd,'rx');

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
                    if C1(i,j)>30
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
                    if C2(i,j)>30
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
    figure(cluster);
    for i=1:length(B)
        if ~isempty(B{i})
            switch i
                case {1}
                    for j=1:length(B{1})
                        plot(XR(B{1}(1,j)),YR(B{1}(1,j)),'ko');
                        hold on;
                    end
                case {2}
                    for j=1:length(B{2})
                        plot(XR(B{2}(1,j)),YR(B{2}(1,j)),'ro');
                        hold on;
                    end
                case {3}
                    for j=1:length(B{3})
                        plot(XR(B{3}(1,j)),YR(B{3}(1,j)),'go');
                        hold on;
                    end    
                 case {4}
                    for j=1:length(B{4})
                        plot(XR(B{4}(1,j)),YR(B{4}(1,j)),'bo');
                        hold on;
                    end 
                 case {5}
                    for j=1:length(B{5})
                        plot(XR(B{5}(1,j)),YR(B{5}(1,j)),'co');
                        hold on;
                    end 
                 case {6}
                    for j=1:length(B{6})
                        plot(XR(B{6}(1,j)),YR(B{6}(1,j)),'yo');
                        hold on;
                    end   
                    case {7}
                    for j=1:length(B{7})
                        plot(XR(B{7}(1,j)),YR(B{7}(1,j)),'mo');
                        hold on;
                    end
                    case {8}
                    for j=1:length(B{8})
                        plot(XR(B{8}(1,j)),YR(B{8}(1,j)),'ko');
                        hold on;
                    end
                case {9}
                    for j=1:length(B{9})
                        plot(XR(B{9}(1,j)),YR(B{9}(1,j)),'ro');
                        hold on;
                    end
                case {10}
                    for j=1:length(B{10})
                        plot(XR(B{10}(1,j)),YR(B{10}(1,j)),'go');
                        hold on;
                    end    
                 case {11}
                    for j=1:length(B{11})
                        plot(XR(B{11}(1,j)),YR(B{11}(1,j)),'bo');
                        hold on;
                    end 
                 case {12}
                    for j=1:length(B{12})
                        plot(XR(B{12}(1,j)),YR(B{12}(1,j)),'co');
                        hold on;
                    end 
                 case {13}
                    for j=1:length(B{13})
                        plot(XR(B{13}(1,j)),YR(B{13}(1,j)),'yo');
                        hold on;
                    end   
                    case {14}
                    for j=1:length(B{14})
                        plot(XR(B{14}(1,j)),YR(B{14}(1,j)),'mo');
                        hold on;
                    end
            end
        end
    end
    plot(S(n+1).xd,S(n+1).yd,'rx');
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
        figure(cluster);        
        plot(S2(B{i}(index{i}(1,1))).xd,S2(B{i}(index{i}(1,1))).yd,'k*');
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

for g=0:1:rmax
    g;
    if(mod(g, round(1/p2) )==0)
        for i=1:1:n
            if S2(i).type ~='L'
                S2(i).G=0;
                S2(i).cl=0;
            end
        end
    end
    
    packets_TO_BS2=0;
    packets_TO_CH2=0;
    
    PACKETS_TO_CH2(g+1)=0;
    PACKETS_TO_BS2(g+1)=0;
    
    for i=1:n
        if (S2(i).E<=0 && S2(i).Dead==0)
            dead2=dead2+1;
            S2(i).Dead=1;
        end
        if S2(i).type~='L' && S2(i).E>0
            S2(i).type='N';
        end
    end
    
    STATISTICS2(g+1).DEAD=dead2;
    DEAD2(g+1)=dead2;
    
    if (dead2==1)
        if(flag_first_dead2==0)
            first_dead2=g;
            flag_first_dead2=1;
        end
    end
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
                        if( temp_rand2 <= (p2/(1-p2*mod(g,round(1/p2)))))
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
                            PACKETS_TO_BS2(g+1) = packets_TO_BS2;
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
                            PACKETS_TO_CH2(g+1) = n - dead2 - cluster2 + 1; %所有的非死亡的普通节点都发送数据包
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
            
            STATISTICS.PACKETS_TO_CH2(g+1)=packets_TO_CH2;
            STATISTICS.PACKETS_TO_BS2(g+1)=packets_TO_BS2;
            countCHs2;
            rcountCHs2 = rcountCHs2 + countCHs2;     
            
            
        end
    end
end


figure(11);
x11=1:1:r;
y11=1:1:r;
x21=1:1:g;
y21=1:1:g;
for i=1:rmax
    x11(i)=i;
    y11(i) = n - STATISTICS1(i).DEAD;
    x21(i)=i;
    y21(i) = n - STATISTICS2(i).DEAD;
end

plot(x11,y11,'r--',x21,y21,'g--');

hold on;
title('图一：存活节点数统计图');
xlabel('工作轮数（/轮）');
ylabel('存活节点数（/个）');
hold off;


figure(12);
x12=1:1:r;
y12=1:1:r;
x22=1:1:r;
y22=1:1:r;

for i=1:rmax
    x12(i)=i;
    y12(i)=STATISTICS.PACKETS_TO_BS1(i+1);
    x22(i)=i;
    y22(i)=STATISTICS.PACKETS_TO_BS2(i+1);
end
plot(x12,y12,'r--',x22,y22,'g--');