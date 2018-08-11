clear;
clc;
N=100;  %�ڵ�ĸ���
sink_x=N;  %sink�ڵ��λ��
sink_y=N;
source_x=0; %Դ�ڵ��λ��
source_y=0;
Eavailable=zeros(1,N);
Einitial=1;
randlocation=N*rand(2,N);  %�������ڵ�
nodelocation=zeros(2,N);
Rt=30;  %���Ͱ뾶
L=4000; % ���ݰ���С Data packet size
Pt=14.88e-003;  %���͹���
Pr=12.50e-003;  %���ܹ���
Ttr=0.4e-003;  %����ÿ�ֽ�����Ҫ��ʱ��
Tre=0.4e-003;  %����ÿ�ֽ�����Ҫ��ʱ��
Packet=5;
Delay=50*(Ttr+Tre);  %�������ݰ�ʱ���ӳ�
T_Delay=zeros(N,N);
%Eelec=ETx=ERx,Radio electronics energy;
Eelec=50*0.000000001;  
ETx=50*0.000000001;
ERx=50*0.000000001;  
%Radio amplifer energy;
Efs=10*0.000000000001;  
Eamp=0.0013*0.000000000001;
%������Դ�ڵ��sink�ڵ����ڵĽڵ�λ�þ���
nodelocation(1,1)=source_x; 
nodelocation(2,1)=source_y;
nodelocation(1,N)=sink_x; 
nodelocation(2,N)=sink_y;
Eavailable(1,1)=100;
Eavailable(1,N)=100;
for i=2:1:N-1
    nodelocation(1,i)=randlocation(1,i-1);
    nodelocation(2,i)=randlocation(2,i-1);
    Eavailable(1,i)=Einitial;
end
%���ڵ�֮��ľ���
D=zeros(N,N);
for i=1:1:N
    for ii=1:1:N
        D(i,ii)=sqrt((nodelocation(1,i)-nodelocation(1,ii))^2+(nodelocation(2,i)-nodelocation(2,ii))^2);
    end
end
%ʱ���ӳپ�����ʽ
for iii=1:1:N
    for iiii=1:1:N
        T_Delay(iii,iiii)=Delay;
    end
end
Hopcount=D(1,N)/Rt;
Ssetpoint=500;
Vavg=zeros(N,N);
T_Delay1=zeros(N,N);
K1=10;
Ttime=0; %������������
Econsume=0; %ƽ����������
Tconsume=0; %ƽ��ʱ���ӳ�
Nconsume=0;  %���������Ľڵ����
Tout=0; %��������
Tconsumeevry=0;
Econsumeevery=0;
M=3000;
for T=1:1:M
    clf;
    for jj=1:1:N
       plot(nodelocation(1,jj),nodelocation(2,jj),'r.','markersize',16)
        hold on
    end
    b=1;
    for a=1:1:N  %��ͬ���裬�γɻ�·
        FS=zeros(N,N);
        Vavg1=zeros(N,N);
        e=zeros(1,N);
        e1=zeros(1,N);
        n=zeros(1,N);
        Uout=zeros(N,N);
        for r=1:1:N
             for u=1:1:N
                 if D(r,u)<=Rt
                     if D(r,u)>0
                         if D(r,N)>D(u,N)
                             FS(r,u)=u;
                             e(r,u)=rand;
                             n(1,r)= n(1,r)+1;
                             e1(1,r)=e1(1,r)+e(r,u);
                             Vavg(r,u)=(D(r,N)-D(u,N))/T_Delay(r,u);
                         end
                     end
                 end
             end
        end
        rr=b;
        for uu=1:1:N
            if FS(rr,uu)~=0
               Vavg(rr,uu)=(D(rr,N)-D(uu,N))/T_Delay(rr,uu); 
               if  Vavg(rr,uu)>Ssetpoint
                   Vavg1(rr,uu)=Vavg(rr,uu);
               end
            end
        end
        a=max(Vavg1(rr,:));
        aa=find(Vavg1(rr,:)==a); 
        if a==0
            Uout(1,rr)=1-K1*(e1(1,rr)/n(1,rr)); 
            m=max(e(rr,:));
            mm=find(e(rr,:)==m);
           if m>=Uout(1,rr)
               Eavailable(1,rr)=Eavailable(1,rr)-Packet*L*(Eelec+Efs*(D(rr,mm)*D(rr,mm)));
               Eavailable(1,mm)=Eavailable(1,mm)-Packet*L*Eelec;
               Econsume=Econsume+Packet*L*(Eelec+Efs*(D(rr,mm)*D(rr,mm)))+Packet*L*Eelec;
              % Tconsume=Tconsume+Delay;
               Nconsume=Nconsume+1;
               XXX=[nodelocation(1,rr),nodelocation(1,mm)];
               YYY=[nodelocation(2,rr),nodelocation(2,mm)];
               %plot(XXX,YYY)
               pp=mm;
           else 
               Tout=Tout+1;
               break;
           end
        else
            if length(aa)<=1
                   Eavailable(1,rr)=Eavailable(1,rr)-Packet*L*(Eelec+Efs*(D(rr,aa)*D(rr,aa)));
                   Eavailable(1,aa)=Eavailable(1,aa)-Packet*L*Eelec;
                   Econsume=Econsume+Packet*L*(Eelec+Efs*(D(rr,aa)*D(rr,aa)))+Packet*L*Eelec;
                   %Tconsume=Tconsume+Delay;
                   Nconsume=Nconsume+1;
                   XXX=[nodelocation(1,rr),nodelocation(1,aa)];
                   YYY=[nodelocation(2,rr),nodelocation(2,aa)];
                   %plot(XXX,YYY)
                   pp=aa;
                else
                    aa=aa(1,1);
                    Eavailable(1,rr)=Eavailable(1,rr)-Packet*L*(Eelec+Efs*(D(rr,aa)*D(rr,aa)));
                    Eavailable(1,aa)=Eavailable(1,aa)-Packet*L*Eelec;
                    Econsume=Econsume+Packet*L*(Eelec+Efs*(D(rr,aa)*D(rr,aa)))+Packet*L*Eelec;
                   % Tconsume=Tconsume+Delay;
                    Nconsume=Nconsume+1;
                    XXX=[nodelocation(1,rr),nodelocation(1,aa)];
                    YYY=[nodelocation(2,rr),nodelocation(2,aa)];
                    %plot(XXX,YYY)
                    pp=aa;
            end
        end
         Tconsume=Tconsume+Delay;
        %pause(0.05)
        if pp==N
             break;
        else   
             b=pp;
        end
    end
     Q=min(Eavailable(1,:));
     QQ=find(Eavailable(1,:)==Q);
     if Q<=0
         break;
     else
       Ttime=Ttime+1;
     end
end
Tconsumeevry=Tconsume/M;
Econsumeevery=Econsume/M;