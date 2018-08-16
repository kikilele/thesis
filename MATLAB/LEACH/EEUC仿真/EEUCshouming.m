clear;
xm=200;
ym=200;% 区域大小
sink.x=0.5*xm;
sink.y=-50;
n=1001;
Eo=0.5;
Rc0=90;
packetLength = 500;%数据包长度 
ctrPacketLength = 500;%控制包长度

%Eelec=Etx=Erx
ETX=50*0.000000001;
ERX=50*0.000000001;
%Transmit Amplifier types
Efs=10*0.000000000001;
Emp=0.0013*0.000000000001;
%Data Aggregation Energy
EDA=5*0.000000001;
do=sqrt(Efs/Emp);
RMAX=4000;


r=0;% 节点i距离BS的距离

figure(1);

for i=1:1:n
%S(i).E=0.5; % 节点初始能量
S(i).xd=rand(1,1)*xm;
XR(i)=S(i).xd;
S(i).yd=rand(1,1)*ym;
YR(i)=S(i).yd;
S(i).r=sqrt((S(i).xd-sink.x)^2+(S(i).yd-sink.y)^2);%节点到基站的距离
  S(i).Rc=(1-0.5*(270-S(i).r)/220)*90;
%plot(S(i).xd,S(i).yd,'o');
  
S(i).type='N';
S(i).E=Eo;
plot(S(i).xd,S(i).yd,'o');
hold on;


end
plot(sink.x,sink.y,'x');

for R=0:1:10
R
% 找出T=0.4时所有的候选节点
B_ID=zeros; % 用来存储候选簇头的ID
p=0;% 统计第一层节点数目
for i=1:1:n
    S(i).T=rand(1,1);
    if S(i).T<0.4
        p=p+1;
        B_ID(p)=i;
        N(p).xd=S(i).xd;
        N(p).yd=S(i).yd;
        N(p).Rc=(1-0.5*(270-S(i).r)/220)*90;  %候选簇头广播竞选消息的半径
        N(p).r=S(i).r;
        N(p).type='B';
        N(p).E=S(i).E;
       
    end
end
% 第一轮选举，节点ID大的当选簇头

%a=0; %用来统计候选簇头的数目
B2_ID=B_ID;
C_ID=zeros;
f=0;

for k=1:1:100

for i=1:1:p
    
   
      if B_ID(i)==max( B_ID)&&B_ID(i)~=0
      % plot(N(i).xd,N(i).yd,'r*','MarkerSize',12); 
     f=f+1;
     C_ID(f)=B_ID(i);
     
       B_ID(i)=0;
       N(i).type='N';
   
          hold on
         for j=1:1:p  % 计算在簇头影响范围内的候选节点
            d=sqrt((N(i).xd-N(j).xd)^2+(N(i).yd-N(j).yd)^2);  
             if d<=N(i).Rc
              N(j).type='N';
             % plot(N(j).xd,N(j).yd,'k+','MarkerSize',12); 
              B_ID(j)=0;
             
            end
         end
      end
  end
if B_ID==0
    break;
end

end


% 减去候选簇头在成簇阶段的能耗，因为开始竞选时，每个候选簇头都发送了竞选消息
 for i=1:1:p
     S(B2_ID(i)).E=S(B2_ID(i)).E-2*( ETX * ctrPacketLength + Efs * ctrPacketLength*( S(B2_ID(i)).Rc*S(B2_ID(i)).Rc));
 end
 
 
 % 当选簇头的节点全网广播自己是簇头的消息
 C_B=zeros;   % 存储簇头到基站的距离
for j=1:1:f
     C_B(j)=S(C_ID(j)).r;
    S( C_ID(j)).E=S( C_ID(j)).E-3*( ETX * ctrPacketLength + Emp * ctrPacketLength* 283*283*283*283);  %可以乘以2，第二个就是簇头把时隙分配结果发送给其他节点,第三个是 簇头间构建路由
end

% 计算全网节点组簇 能耗,数据发送能耗
for i=1:1:n
    D=zeros;  
    for j=1:1:f
        d_CH=sqrt((S(i).xd-S(C_ID(j)).xd)^2+(S(i).yd-S(C_ID(j)).yd)^2); %计算节点和 各个簇头的距离 并存入D中
        D(j)=d_CH;
    end
    d1=min(D);
    if d1~=0&&d1<do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Efs * packetLength*d1*d1);
    elseif d1>do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Emp * packetLength*d1*d1*d1*d1);
    end
    u=find(D==min(D));  % 找出最小值在数组中的位置
    S( C_ID(u)).E= S( C_ID(u)).E-ERX* (ctrPacketLength+packetLength)-EDA*packetLength; % 簇头节点减去接受和融合簇内节点信息的能耗
end

% 计算簇头发送以及 转发数据能耗
C_B2=sort(C_B);  %把簇头到 基站的距离 从小到大排列
  C_ID2=zeros;   %把ID号 也是按照r从小到大 排列
   for j=1:1:f
       for k=1:1:f
           if C_B2(j)==S(C_ID(k)).r
              C_ID2(j)= C_ID(k);
           end
       end
   end
% 计算 最远簇头与其他簇头的距离
C_C=cell(f,1); % 存入其他簇头距离第j个簇头的距离
 for j=f:-1:2
      
    for i=1:1:f
      if S(C_ID2(j)).r>S(C_ID2(i)).r
        d=sqrt((S(C_ID2(j)).xd-S(C_ID2(i)).xd)^2+(S(C_ID2(j)).yd-S(C_ID2(i)).yd)^2); 
        C_C{j}(i)=d;
      end
    end   
 end    %完毕
     
 for j=f:-1:1
    if S(C_ID2(j)).r>140
     e=ceil(min(C_C{j})/S(C_ID2(j)).Rc);
    if e*S(C_ID2(j)).Rc<do
         S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Efs * packetLength*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc));
    else
         S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Emp * packetLength*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc));
    end
      g=find(C_C{j}==min(C_C{j}));
      e2=ceil(min(C_C{g})/S(C_ID2(g)).Rc);
      
      if e*S(C_ID2(g)).Rc<do
          S(C_ID2(g)).E= S(C_ID2(g)).E-( ETX * packetLength + Efs * packetLength*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc))- ERX * packetLength;
      else
          S(C_ID2(g)).E= S(C_ID2(g)).E-( ETX * packetLength + Emp * packetLength*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc))- ERX * packetLength;
      end
      
      g2=find(C_C{g}==min(C_C{g}));
     if g2>1
      if S(C_ID2(g2)).r<do
          S(C_ID2(g2)).E= S(C_ID2(g2)).E-( ETX * packetLength + Efs * packetLength*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc))- ERX * packetLength;
           else
          S(C_ID2(g2)).E= S(C_ID2(g2)).E-( ETX * packetLength + Emp * packetLength*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc))- ERX * packetLength;
      end
     end
    else
        if S(C_ID2(j)).r<do
            S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Efs * packetLength*S(C_ID2(j)).r*S(C_ID2(j)).r); 
        else
           S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Emp * packetLength*S(C_ID2(j)).r*S(C_ID2(j)).r*S(C_ID2(j)).r*S(C_ID2(j)).r);
        end
    end
 end
end


DEAD=zeros;
LIVE=zeros;
for R=10:1:RMAX
    R
    dead=0;
B_ID=zeros; % 用来存储候选簇头的ID
E=zeros;
p=0;% 统计第一层节点数目
for i=1:1:n
  if S(i).E>0
    S(i).T=rand(1,1);
    if S(i).T<0.4
        p=p+1;
        B_ID(p)=i;
        E(p)=S(i).E;
        N(p).xd=S(i).xd;
        N(p).yd=S(i).yd;
        N(p).Rc=(1-0.5*(270-S(i).r)/220)*90;  %候选簇头广播竞选消息的半径
        N(p).r=S(i).r;
        N(p).type='B';
        N(p).E=S(i).E;
       
    end
  end
end

if p==0
    break;
end


%a=0; %用来统计候选簇头的数目
B2_ID=B_ID;
C_ID=zeros;
f=0;

for k=1:1:100


    
      z=find(E==max(E));
      
      % plot(N(z).xd,N(z).yd,'r*','MarkerSize',12); 
     f=f+1;
     C_ID(f)=B_ID(z);
     E(z)=0;
       B_ID(z)=0;
       N(z).type='N';
   
          hold on
         for j=1:1:p  % 计算在簇头影响范围内的候选节点
            d=sqrt((N(z).xd-N(j).xd)^2+(N(z).yd-N(j).yd)^2);  
             if d<=N(z).Rc
              N(j).type='N';
             % plot(N(j).xd,N(j).yd,'k+','MarkerSize',12); 
              B_ID(j)=0;
               E(j)=0;
             
            end
         end
     
     if B_ID==0
       break;
     end

end


% 减去候选簇头在成簇阶段的能耗，因为开始竞选时，每个候选簇头都发送了竞选消息
 for i=1:1:p
     S(B2_ID(i)).E=S(B2_ID(i)).E-2*( ETX * ctrPacketLength + Efs * ctrPacketLength*( S(B2_ID(i)).Rc*S(B2_ID(i)).Rc));
 end
 
 
 % 当选簇头的节点全网广播自己是簇头的消息
 C_B=zeros;   % 存储簇头到基站的距离
for j=1:1:f
     C_B(j)=S(C_ID(j)).r;
    S( C_ID(j)).E=S( C_ID(j)).E-3*( ETX * ctrPacketLength + Emp * ctrPacketLength* 283*283*283*283);  %可以乘以2，第二个就是簇头把时隙分配结果发送给其他节点,第三个是 簇头间构建路由
end

% 计算全网节点组簇 能耗,数据发送能耗
for i=1:1:n
    if S(i).E>0
    D=zeros;  
    for j=1:1:f
        d_CH=sqrt((S(i).xd-S(C_ID(j)).xd)^2+(S(i).yd-S(C_ID(j)).yd)^2); %计算节点和 各个簇头的距离 并存入D中
        D(j)=d_CH;
    end
    d1=min(D);
    if d1~=0&&d1<do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Efs * packetLength*d1*d1);
    elseif d1>do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Emp * packetLength*d1*d1*d1*d1);
    end
    u=find(D==min(D));  % 找出最小值在数组中的位置
    S( C_ID(u)).E= S( C_ID(u)).E-ERX* (ctrPacketLength+packetLength)-EDA*packetLength; % 簇头节点减去接受和融合簇内节点信息的能耗
    end
end

% 计算簇头发送以及 转发数据能耗
C_B2=sort(C_B);  %把簇头到 基站的距离 从小到大排列
  C_ID2=zeros;   %把ID号 也是按照r从小到大 排列
   for j=1:1:f
       for k=1:1:f
           if C_B2(j)==S(C_ID(k)).r
              C_ID2(j)= C_ID(k);
           end
       end
   end
% 计算 最远簇头与其他簇头的距离
C_C=cell(f,1); % 存入其他簇头距离第j个簇头的距离
 for j=f:-1:2
      
    for i=1:1:f
      if S(C_ID2(j)).r>S(C_ID2(i)).r
        d=sqrt((S(C_ID2(j)).xd-S(C_ID2(i)).xd)^2+(S(C_ID2(j)).yd-S(C_ID2(i)).yd)^2); 
        C_C{j}(i)=d;
      end
    end   
 end    %完毕
     
 for j=f:-1:1
    if S(C_ID2(j)).r>140
     e=ceil(min(C_C{j})/S(C_ID2(j)).Rc);
    if e*S(C_ID2(j)).Rc<do
         S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Efs * packetLength*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc));
    else
         S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Emp * packetLength*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc)*(e*S(C_ID2(j)).Rc));
    end
    if g>1
    g=find(C_C{j}==min(C_C{j}));
      e2=ceil(min(C_C{g})/S(C_ID2(g)).Rc);
      
      if e*S(C_ID2(g)).Rc<do
          S(C_ID2(g)).E= S(C_ID2(g)).E-( ETX * packetLength + Efs * packetLength*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc))- ERX * packetLength;
      else
          S(C_ID2(g)).E= S(C_ID2(g)).E-( ETX * packetLength + Emp * packetLength*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc)*(e*S(C_ID2(g)).Rc))- ERX * packetLength;
      end
    end
      
      g2=find(C_C{g}==min(C_C{g}));
     if g2>1
      if S(C_ID2(g2)).r<do
          S(C_ID2(g2)).E= S(C_ID2(g2)).E-( ETX * packetLength + Efs * packetLength*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc))- ERX * packetLength;
           else
          S(C_ID2(g2)).E= S(C_ID2(g2)).E-( ETX * packetLength + Emp * packetLength*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc)*(e*S(C_ID2(g2)).Rc))- ERX * packetLength;
      end
     end
    else
        if S(C_ID2(j)).r<do
            S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Efs * packetLength*S(C_ID2(j)).r*S(C_ID2(j)).r); 
        else
           S(C_ID2(j)).E=S(C_ID2(j)).E-( ETX * packetLength + Emp * packetLength*S(C_ID2(j)).r*S(C_ID2(j)).r*S(C_ID2(j)).r*S(C_ID2(j)).r);
        end
    end
 end
 
 for i=1:1:n
    if S(i).E<=0
        dead=dead+1;
        
    end
 end
 DEAD(R)=dead;
        LIVE(R-9)=1001-dead;
 if dead>995
     break;
 end
end



%axis([-50,250,-50,250])
%axis equal   %纵横坐标采用等长刻度     axis 设置做标轴的最小最大值  用法  axis([xmin,xmax,ymin,ymax])
%hold on
%x=[0 200 200 0 0];
%y=[0 0 200 200 0];
%plot(x,y,'--')
%hold on;
%xlabel('X/m','fontname','宋体','fontsize',12);
%ylabel('Y/m','fontname','宋体','fontsize',12);
%box off
    