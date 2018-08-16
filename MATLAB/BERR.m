clear;
clc;
nodenumber=100;
minimumreliablebound=0.9;
length=nodenumber;
width=nodenumber;
sink_x=nodenumber;  %sink�ڵ��λ��
sink_y=nodenumber;
source_x=0; %Դ�ڵ��λ��
source_y=0;
sinkposition=[length;width];
energyinitialvalue=1;
circle=10;
reliabletransmissiondistance=30;
packetnumber=21;  %���ݰ�����
packetlength=4000; % ���ݰ���С Data packet size
% delay��ȷ��
Ttr=0.4e-003;  %����ÿ�ֽ�����Ҫ��ʱ��
Tre=0.4e-003;  %����ÿ�ֽ�����Ҫ��ʱ��
Delay=50*(Ttr+Tre);  %�������ݰ�ʱ���ӳ�
%Radio amplifer energy;
Eelec=50*0.000000001;
Efs=10*0.000000000001;  
Eamp=0.0013*0.000000000001;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ڵ�λ����Ϣ
%������Դ�ڵ��sink�ڵ����ڵĽڵ�λ�þ���
nodeposition(1,1)=source_x; 
nodeposition(2,1)=source_y;
nodeposition(1,nodenumber)=sink_x; 
nodeposition(2,nodenumber)=sink_y;
randlocation=nodenumber*rand(2,nodenumber);  %�������ڵ�
nodeenergyinitial(1,1)=100;
nodeenergyinitial(1,nodenumber)=100;
for i=2:1:nodenumber-1
    nodeposition(1,i)=randlocation(1,i-1);
    nodeposition(2,i)=randlocation(2,i-1);
    plot(nodeposition(1,i),nodeposition(2,i),'r.','markersize',16);
    nodeenergyinitial(1,i)=energyinitialvalue;
    hold on;
end
plot(nodeposition(1,1),nodeposition(2,1),'r^-','markersize',16);
plot(sinkposition(1,1),sinkposition(2,1),'r^-','markersize',16);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ѱ��һ���ھӽڵ��б�
%���������֮��ľ���
onehoplist=zeros(nodenumber,nodenumber); %zeros(100,100)
for ii=1:nodenumber
   for j=1:1:nodenumber
       D(j,ii)=sqrt((nodeposition(1,ii)-nodeposition(1,j))^2+(nodeposition(2,ii)-nodeposition(2,j))^2); %ŷʽ����
    end
end
for i=1:nodenumber
    for j=1:nodenumber
        if  i~=j
            if D(i,j)>0 && D(i,j)<=reliabletransmissiondistance   %30
                onehoplist(j,i)=j;   
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%sink�ڵ��һ���ھӽڵ��б�
hopcount=zeros(1,nodenumber);
hoplist=zeros(nodenumber,nodenumber);   %sink�ڵ��n����Ϣ��
for sink=nodenumber
    ki=1;
    for jj=1:nodenumber
       if jj~=sink
          if D(sink,jj)<=reliabletransmissiondistance
             hoplist(1,jj)=jj;    %sink�ڵ��һ����Ϣ��
             hopcount(1,jj)=ki;
             text(nodeposition(1,jj)+1,nodeposition(2,jj)+1,'HC1');
          end
        end
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ѱ�������ھӽڵ��б�
  for i=nodenumber
      ki=2;
      for j=1:nodenumber
        for k=1:nodenumber
            if onehoplist(j,i)~=0 && onehoplist(k,j)~=0 && onehoplist(k,j)~=i && hopcount(1,onehoplist(k,j))==0
               hopcount(1,k)=ki;     
               hoplist(ki,k)=k;
               text(nodeposition(1,k)+1,nodeposition(2,k)+1,'HC2');
            end
        end
    end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ѱ�������ھӽڵ��б�
 for i=nodenumber
      ki=3;
      for j=1:nodenumber
        for k=1:nodenumber
            if hopcount(1,j)==2 && onehoplist(k,j)~=0 && hopcount(1,k)==0
               hopcount(1,k)=ki;     
               hoplist(ki,k)=k;
               text(nodeposition(1,k)+1,nodeposition(2,k)+1,'HC3');
            end
        end
    end
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ѱ�������ھӽڵ��б�
 for i=nodenumber
      ki=4;
      for j=1:nodenumber
        for k=1:nodenumber
            if hopcount(1,j)==3 && onehoplist(k,j)~=0 && hopcount(1,k)==0
               hopcount(1,k)=ki;     
               hoplist(ki,k)=k;
               text(nodeposition(1,k)+1,nodeposition(2,k)+1,'HC4');
            end
        end
    end
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ѱ�������ھӽڵ��б�
 for i=nodenumber
      ki=5;
      for j=1:nodenumber
        for k=1:nodenumber
            if hopcount(1,j)==4 && onehoplist(k,j)~=0 && hopcount(1,k)==0
               hopcount(1,k)=ki;     
               hoplist(ki,k)=k;
               text(nodeposition(1,k)+1,nodeposition(2,k)+1,'HC5');
            end
        end
    end
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Ѱ�������ھӽڵ��б�
 for i=nodenumber
      ki=6;
      for j=1:nodenumber
        for k=1:nodenumber
            if hopcount(1,j)==5 && onehoplist(k,j)~=0 && hopcount(1,k)==0
               hopcount(1,k)=ki;     
               hoplist(ki,k)=k;
               text(nodeposition(1,k)+1,nodeposition(2,k)+1,'HC6');
            end
        end
    end
 end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%���������ֵ
HCmax=max(hopcount(1,:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%����·�����ʣ��ھ��б�Energyefficiency
Renumber=zeros(nodenumber,nodenumber);
p=zeros(nodenumber,nodenumber);
F=zeros(nodenumber,nodenumber);
for q1=1:nodenumber
    for q2=1:nodenumber
        if onehoplist(q2,q1)~=0 && D(q1,nodenumber)>=D(q2,nodenumber)
           p(q1,q2)=rand;   %�������������·���ݰ�������
           F(q1,q2)=q2;      %�ھӽڵ��
           Renumber(q1,q2)=1/(1-p(q1,q2));  %�ɹ��ɿ��ķ���һ�����ݰ���ϣ���ĸ���
        end
    end
end
a=0.1;
Parthdelay=0;     %��ѡ·�������ӳ�
Econsumepathtotal=0;   %��ѡ·��������������
Nconsume=0;   %���õ��ڵ�ĸ���
Ttime=0;
Parthdelayevery=0;
Econsumepathtotaevery=0;
M=100;    %����
for T=1:1:M
    clf;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ڵ�λ����Ϣ
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%������Դ�ڵ��sink�ڵ����ڵĽڵ�λ�þ���
    for i=2:1:nodenumber-1
        %nodeposition(1,i)=randlocation(1,i-1);
        %nodeposition(2,i)=randlocation(2,i-1);
        plot(nodeposition(1,i),nodeposition(2,i),'r.','markersize',16);
        hold on;
    end
    plot(nodeposition(1,1),nodeposition(2,1),'r^-','markersize',16);
    plot(sinkposition(1,1),sinkposition(2,1),'r^-','markersize',16);
   q5=1;
   for Y=1:1:nodenumber 
       ECopt=zeros(1,nodenumber);
       clear ABility;
       EC=zeros(nodenumber,nodenumber);
       Econsume=zeros(nodenumber,nodenumber);
       Ptradeoff=zeros(nodenumber,nodenumber);
       PEC=zeros(nodenumber,nodenumber);
       for q3=1:nodenumber-1
         for q4=1:nodenumber
             if F(q3,q4)~=0
                 if q4~=nodenumber
                    Econsume(q3,q4)=packetnumber*packetlength*(Eelec+Efs*(D(q3,q4)*D(q3,q4)));
                    EC(q3,q4)=Econsume(q3,q4)*Renumber(q3,q4)/nodeenergyinitial(1,q3);
                    %  ABility(q3,q4)=1/EC(q3,q4);
                    %ECtotal(1,q3)=ECtotal(1,q3)+EC(q3,q4);
                 else
                    Econsume(q3,q4)=packetnumber*packetlength*(Eelec+Efs*(D(q3,q4)*D(q3,q4)));
                    EC(q3,q4)=Econsume(q3,q4)/nodeenergyinitial(1,q3);
                    %ABility(q3,q4)=1/EC(q3,q4);
                   % ECtotal(1,q3)=ECtotal(1,q3)+EC(q3,q4);
                 end
             end
         end
        % w=max(ABility(q3,:));
         %ww=find(ABility(q3,:)==w);
         %ECopt(1,q3)=EC(q3,ww);
          w=find(EC(q3,:)==0);
         [ww www]=size(w);
         for qq=1:1:www
            EC(q3,w(qq))=1000;
         end
         y=min(EC(q3,:));
         yy=find(EC(q3,:)==y);
         ECopt(1,q3)=EC(q3,yy);
      end
       for q6=1:nodenumber
           if F(q5,q6)~=0
              PEC(q5,q6)=EC(q5,q6)+ECopt(1,q6);
              Ptradeoff(q5,q6)=a*PEC(q5,q6)+(1-a)*hopcount(1,q6)/HCmax;  
           end
       end
       h=find(Ptradeoff(q5,:)==0);
       [hh hhh]=size(h);
       for b=1:1:hhh
           Ptradeoff(q5,h(b))=1000;
       end
       u=min(Ptradeoff(q5,:));
       uu=find(Ptradeoff(q5,:)==u);
       nodeenergyinitial(1,q5)=nodeenergyinitial(1,q5)-packetnumber*packetlength*(Eelec+Efs*(D(q3,q4)*D(q3,q4)));  %����ڵ��ʣ����������
       nodeenergyinitial(1,uu)=nodeenergyinitial(1,uu)-packetnumber*packetlength*Eelec;  %����ڵ��ʣ����������
       Econsumepathtotal=Econsumepathtotal+Econsume(q5,uu)+packetnumber*packetlength*Eelec;   %��ѡ·��������������
       Parthdelay=Parthdelay+Delay;   %��ѡ·�������ӳ�
       Nconsume=Nconsume+1;    %���õ��Ľڵ����
       XXX=[nodeposition(1,q5),nodeposition(1,uu)];
       YYY=[nodeposition(2,q5),nodeposition(2,uu)];
       %plot(XXX,YYY)
       q5=uu;
       %pause(0.05)
       if  q5==nodenumber
           break;
       end
   end
   Q=min(nodeenergyinitial(1,:));
   QQ=find(nodeenergyinitial(1,:)==Q);
   if Q<=0
      break;
   else
         Ttime=Ttime+1;
   end
end
Parthdelayevery=Parthdelay/Ttime;
Econsumepathtotaevery=Econsumepathtotal/Ttime;