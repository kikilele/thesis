%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  DV-Hop�㷨  ~~~~~~~~~~~~~~~~~~~~~~~~
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% BorderLength-----����������ı߳�����λ��m
% NodeAmount-------����ڵ�ĸ���
% BeaconAmount---�ű�ڵ���
% Sxy--------------���ڴ洢�ڵ����ţ������꣬������ľ���
%Beacon----------�ű�ڵ��������;BeaconAmount*BeaconAmount
%UN-------------δ֪�ڵ��������;2*UNAmount
% Distance------δ֪�ڵ㵽�ű�ڵ�������;2*BeaconAmount
%h---------------�ڵ���ʼ��������
%X---------------�ڵ���������ʼ����,X=[x,y]'
% R------------------�ڵ��ͨ�ž��룬һ��Ϊ10-100m

clear,close all;
BorderLength=100;
NodeAmount=100;
BeaconAmount=8;
UNAmount=NodeAmount-BeaconAmount;
R=50;
% D=zeros(NodeAmount,NodeAmount);%δ֪�ڵ絽�ű�ڵ�����ʼ����BeaconAmount��NodeAmount��
h=zeros(NodeAmount,NodeAmount);%��ʼ����Ϊ0��BeaconAmount��NodeAmount��
X=zeros(2,UNAmount);%�ڵ���������ʼ����

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~�������������ڲ������ȷֲ����������~~~~~~~~~~~~~~~~~~~~
C=BorderLength.*rand(2,NodeAmount);
%���߼��ŵĽڵ�����
Sxy=[[1:NodeAmount];C];
Beacon=[Sxy(2,1:BeaconAmount);Sxy(3,1:BeaconAmount)];%�ű�ڵ�����
UN=[Sxy(2,(BeaconAmount+1):NodeAmount);Sxy(3,(BeaconAmount+1):NodeAmount)];%δ֪�ڵ�����
%�����ڵ�ֲ�ͼ
plot(Sxy(2,1:BeaconAmount),Sxy(3,1:BeaconAmount),'r*',Sxy(2,(BeaconAmount+1):NodeAmount),Sxy(3,(BeaconAmount+1):NodeAmount),'k.')
xlim([0,BorderLength]);
ylim([0,BorderLength]);
title('* ��ɫ�ű�ڵ� . ��ɫδ֪�ڵ�')
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~��ʼ���ڵ����롢��������~~~~~~~~~~~~~~~~~~~~~~
for i=1:NodeAmount
    for j=1:NodeAmount
        Dall(i,j)=((Sxy(2,i)-Sxy(2,j))^2+(Sxy(3,i)-Sxy(3,j))^2)^0.5;%���нڵ���໥����
        if (Dall(i,j)<=R)&(Dall(i,j)>0)
            h(i,j)=1;%��ʼ��������
        elseif i==j
            h(i,j)=0;
        else h(i,j)=inf;
        end
    end
end
%~~~~~~~~~~~~~~~~~~~~~~~~~���·���㷨����ڵ������~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for k=1:NodeAmount
    for i=1:NodeAmount
        for j=1:NodeAmount
            if h(i,k)+h(k,j)<h(i,j)%min(h(i,j),h(i,k)+h(k,j))
                h(i,j)=h(i,k)+h(k,j);
            end
        end
    end
end
h
%~~~~~~~~~~~~~~~~~~~~~~~~~��ÿ���ű�ڵ��У��ֵ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
h1=h(1:BeaconAmount,1:BeaconAmount); 
D1=Dall(1:BeaconAmount,1:BeaconAmount);
for i=1:BeaconAmount
    dhop(i,1)=sum(D1(i,:))/sum(h1(i,:));%ÿ���ű�ڵ��ƽ��ÿ������
end
D2=Dall(1:BeaconAmount,(BeaconAmount+1):NodeAmount);%BeaconAmount��UNAmount��
for i=1:BeaconAmount
    for j=1:UNAmount
        if min(D2(:,j))==D2(i,j)
            Dhop(1,j)=D2(i,j);%δ֪�ڵ��������ű���У��ֵ
        end
    end
end
Dhop
%~~~~~~~~~~~~~~~~~~~~~~~~~~~���������ƾ���~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
hop1=h(1:BeaconAmount,(BeaconAmount+1):NodeAmount)%δ֪�ڵ㵽�ű�������BeaconAmount��UNAmount��
for i=1:UNAmount
    hop=Dhop(1,i);%hopΪ������ű��õ�У��ֵ
    Distance(:,i)=hop*hop1(:,i);%%Beacon��UN�У�
end
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~��С���˷���δ֪������~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d=Distance;
for i=1:2
    for j=1:(BeaconAmount-1)
      a(i,j)=Beacon(i,j)-Beacon(i,BeaconAmount);
    end
end
A=-2*(a');
% d=d1';
 for m=1:UNAmount 
     for i=1:(BeaconAmount-1)
         B(i,1)=d(i,m)^2-d(BeaconAmount,m)^2-Beacon(1,i)^2+Beacon(1,BeaconAmount)^2-Beacon(2,i)^2+Beacon(2,BeaconAmount)^2;
     end
           X1=inv(A'*A)*A'*B;
           X(1,m)=X1(1,1);
           X(2,m)=X1(2,1);
 end
 UN
 X
 for i=1:UNAmount
     error(1,i)=(((X(1,i)-UN(1,i))^2+(X(2,i)-UN(2,i))^2)^0.5);
 end
 figure;plot(error,'-o')
 title('ÿ��δ֪�ڵ�����')
 error=sum(error)/UNAmount
 Accuracy=error/R