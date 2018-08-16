clear;
xm=200;
ym=200;% �����С
sink.x=0.5*xm;
sink.y=-50;
n=1001;
Eo=0.5;
Rc0=90;
packetLength = 500;%���ݰ����� 
ctrPacketLength = 500;%���ư�����

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


r=0;% �ڵ�i����BS�ľ���

figure(1);

for i=1:1:n
%S(i).E=0.5; % �ڵ��ʼ����
S(i).xd=rand(1,1)*xm;
XR(i)=S(i).xd;
S(i).yd=rand(1,1)*ym;
YR(i)=S(i).yd;
S(i).r=sqrt((S(i).xd-sink.x)^2+(S(i).yd-sink.y)^2);%�ڵ㵽��վ�ľ���
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
% �ҳ�T=0.4ʱ���еĺ�ѡ�ڵ�
B_ID=zeros; % �����洢��ѡ��ͷ��ID
p=0;% ͳ�Ƶ�һ��ڵ���Ŀ
for i=1:1:n
    S(i).T=rand(1,1);
    if S(i).T<0.4
        p=p+1;
        B_ID(p)=i;
        N(p).xd=S(i).xd;
        N(p).yd=S(i).yd;
        N(p).Rc=(1-0.5*(270-S(i).r)/220)*90;  %��ѡ��ͷ�㲥��ѡ��Ϣ�İ뾶
        N(p).r=S(i).r;
        N(p).type='B';
        N(p).E=S(i).E;
       
    end
end
% ��һ��ѡ�٣��ڵ�ID��ĵ�ѡ��ͷ

%a=0; %����ͳ�ƺ�ѡ��ͷ����Ŀ
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
         for j=1:1:p  % �����ڴ�ͷӰ�췶Χ�ڵĺ�ѡ�ڵ�
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


% ��ȥ��ѡ��ͷ�ڳɴؽ׶ε��ܺģ���Ϊ��ʼ��ѡʱ��ÿ����ѡ��ͷ�������˾�ѡ��Ϣ
 for i=1:1:p
     S(B2_ID(i)).E=S(B2_ID(i)).E-2*( ETX * ctrPacketLength + Efs * ctrPacketLength*( S(B2_ID(i)).Rc*S(B2_ID(i)).Rc));
 end
 
 
 % ��ѡ��ͷ�Ľڵ�ȫ���㲥�Լ��Ǵ�ͷ����Ϣ
 C_B=zeros;   % �洢��ͷ����վ�ľ���
for j=1:1:f
     C_B(j)=S(C_ID(j)).r;
    S( C_ID(j)).E=S( C_ID(j)).E-3*( ETX * ctrPacketLength + Emp * ctrPacketLength* 283*283*283*283);  %���Գ���2���ڶ������Ǵ�ͷ��ʱ϶���������͸������ڵ�,�������� ��ͷ�乹��·��
end

% ����ȫ���ڵ���� �ܺ�,���ݷ����ܺ�
for i=1:1:n
    D=zeros;  
    for j=1:1:f
        d_CH=sqrt((S(i).xd-S(C_ID(j)).xd)^2+(S(i).yd-S(C_ID(j)).yd)^2); %����ڵ�� ������ͷ�ľ��� ������D��
        D(j)=d_CH;
    end
    d1=min(D);
    if d1~=0&&d1<do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Efs * packetLength*d1*d1);
    elseif d1>do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Emp * packetLength*d1*d1*d1*d1);
    end
    u=find(D==min(D));  % �ҳ���Сֵ�������е�λ��
    S( C_ID(u)).E= S( C_ID(u)).E-ERX* (ctrPacketLength+packetLength)-EDA*packetLength; % ��ͷ�ڵ��ȥ���ܺ��ںϴ��ڽڵ���Ϣ���ܺ�
end

% �����ͷ�����Լ� ת�������ܺ�
C_B2=sort(C_B);  %�Ѵ�ͷ�� ��վ�ľ��� ��С��������
  C_ID2=zeros;   %��ID�� Ҳ�ǰ���r��С���� ����
   for j=1:1:f
       for k=1:1:f
           if C_B2(j)==S(C_ID(k)).r
              C_ID2(j)= C_ID(k);
           end
       end
   end
% ���� ��Զ��ͷ��������ͷ�ľ���
C_C=cell(f,1); % ����������ͷ�����j����ͷ�ľ���
 for j=f:-1:2
      
    for i=1:1:f
      if S(C_ID2(j)).r>S(C_ID2(i)).r
        d=sqrt((S(C_ID2(j)).xd-S(C_ID2(i)).xd)^2+(S(C_ID2(j)).yd-S(C_ID2(i)).yd)^2); 
        C_C{j}(i)=d;
      end
    end   
 end    %���
     
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
B_ID=zeros; % �����洢��ѡ��ͷ��ID
E=zeros;
p=0;% ͳ�Ƶ�һ��ڵ���Ŀ
for i=1:1:n
  if S(i).E>0
    S(i).T=rand(1,1);
    if S(i).T<0.4
        p=p+1;
        B_ID(p)=i;
        E(p)=S(i).E;
        N(p).xd=S(i).xd;
        N(p).yd=S(i).yd;
        N(p).Rc=(1-0.5*(270-S(i).r)/220)*90;  %��ѡ��ͷ�㲥��ѡ��Ϣ�İ뾶
        N(p).r=S(i).r;
        N(p).type='B';
        N(p).E=S(i).E;
       
    end
  end
end

if p==0
    break;
end


%a=0; %����ͳ�ƺ�ѡ��ͷ����Ŀ
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
         for j=1:1:p  % �����ڴ�ͷӰ�췶Χ�ڵĺ�ѡ�ڵ�
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


% ��ȥ��ѡ��ͷ�ڳɴؽ׶ε��ܺģ���Ϊ��ʼ��ѡʱ��ÿ����ѡ��ͷ�������˾�ѡ��Ϣ
 for i=1:1:p
     S(B2_ID(i)).E=S(B2_ID(i)).E-2*( ETX * ctrPacketLength + Efs * ctrPacketLength*( S(B2_ID(i)).Rc*S(B2_ID(i)).Rc));
 end
 
 
 % ��ѡ��ͷ�Ľڵ�ȫ���㲥�Լ��Ǵ�ͷ����Ϣ
 C_B=zeros;   % �洢��ͷ����վ�ľ���
for j=1:1:f
     C_B(j)=S(C_ID(j)).r;
    S( C_ID(j)).E=S( C_ID(j)).E-3*( ETX * ctrPacketLength + Emp * ctrPacketLength* 283*283*283*283);  %���Գ���2���ڶ������Ǵ�ͷ��ʱ϶���������͸������ڵ�,�������� ��ͷ�乹��·��
end

% ����ȫ���ڵ���� �ܺ�,���ݷ����ܺ�
for i=1:1:n
    if S(i).E>0
    D=zeros;  
    for j=1:1:f
        d_CH=sqrt((S(i).xd-S(C_ID(j)).xd)^2+(S(i).yd-S(C_ID(j)).yd)^2); %����ڵ�� ������ͷ�ľ��� ������D��
        D(j)=d_CH;
    end
    d1=min(D);
    if d1~=0&&d1<do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Efs * packetLength*d1*d1);
    elseif d1>do
        S(i).E=S(i).E-2*ERX* ctrPacketLength-( ETX * packetLength + Emp * packetLength*d1*d1*d1*d1);
    end
    u=find(D==min(D));  % �ҳ���Сֵ�������е�λ��
    S( C_ID(u)).E= S( C_ID(u)).E-ERX* (ctrPacketLength+packetLength)-EDA*packetLength; % ��ͷ�ڵ��ȥ���ܺ��ںϴ��ڽڵ���Ϣ���ܺ�
    end
end

% �����ͷ�����Լ� ת�������ܺ�
C_B2=sort(C_B);  %�Ѵ�ͷ�� ��վ�ľ��� ��С��������
  C_ID2=zeros;   %��ID�� Ҳ�ǰ���r��С���� ����
   for j=1:1:f
       for k=1:1:f
           if C_B2(j)==S(C_ID(k)).r
              C_ID2(j)= C_ID(k);
           end
       end
   end
% ���� ��Զ��ͷ��������ͷ�ľ���
C_C=cell(f,1); % ����������ͷ�����j����ͷ�ľ���
 for j=f:-1:2
      
    for i=1:1:f
      if S(C_ID2(j)).r>S(C_ID2(i)).r
        d=sqrt((S(C_ID2(j)).xd-S(C_ID2(i)).xd)^2+(S(C_ID2(j)).yd-S(C_ID2(i)).yd)^2); 
        C_C{j}(i)=d;
      end
    end   
 end    %���
     
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
%axis equal   %�ݺ�������õȳ��̶�     axis �������������С���ֵ  �÷�  axis([xmin,xmax,ymin,ymax])
%hold on
%x=[0 200 200 0 0];
%y=[0 0 200 200 0];
%plot(x,y,'--')
%hold on;
%xlabel('X/m','fontname','����','fontsize',12);
%ylabel('Y/m','fontname','����','fontsize',12);
%box off
    