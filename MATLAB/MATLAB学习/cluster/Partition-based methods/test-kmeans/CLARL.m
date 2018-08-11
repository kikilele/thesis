clc;  
clear;  
  
load Data3.mat;  
  
k=3; %�����������Ŀ  
  
time=5;%timeΪ�����Ĵ���  
number=30;%numberΪ����������  
for T=1:time  
    ClomStaticSample=zeros(1,number);  
    ClomStaticSample=randsample(ClomStatic,number);   %ClomStaticSample������������  
                                                      %����������������ʹ��kmediod�㷨���о���  
                                                        
    %����������������������������  
    p=randperm(number);  
    Temp=p(1:k);  
    Center=zeros(1,k);  
    for j=1:k  
        Center(j)=ClomStaticSample(Temp(j));  
    end  
    [ClomStaticSample]=sort(ClomStaticSample);  
      
    TempDistance=zeros(number,3);           %�ݴ��ֵ  
      
     while 1  
        Circulm=1;                          %ѭ������  
  
        p1=1;  
        p2=1;  
        p3=1;  
  
        if(Circulm~=1)  
            clear Group1 Group2 Group3;     
        end  
        for i=1:number  
            for j=1:3  
                TempDistance(i,j)=abs(ClomStaticSample(i)-Center(j));  
            end  
            [RowMin, RowIndex]=min(TempDistance(i,:));  
            if(RowIndex(1)==1)  
                Group1(p1)=ClomStaticSample(i);  
                p1=p1+1;  
            elseif(RowIndex(1)==2)  
                Group2(p2)=ClomStaticSample(i);  
                p2=p2+1;  
            elseif(RowIndex(1)==3)  
                Group3(p3)=ClomStaticSample(i);  
                p3=p3+1;  
            end  
        end  
  
            len1=length(Group1);  
            len2=length(Group2);  
            len3=length(Group3);  
              
  
                  %�ֱ����ÿ�����г��������ĵĵ㵽�������е�ľ����E,E��СʱΪ�����µľ�������.  
                  E=zeros(1,len1-1);  
                  q1=1;  
                  for j=1:len1  
                      for i=1:number  
                        if(Group1(j)~=Center(1)&&i~=j)  
                            E(q1)=floor(abs(Group1(j)-ClomStaticSample(i)));  
                            q1=q1+1;  
                        end  
                      end  
                  end  
                  NewCenter(1)=min(E);  
  
                 E=zeros(1,len2-1);  
                  q2=1;  
                  for j=1:len2  
                      for i=1:number  
                        if(Group2(j)~=Center(2)&&i~=j)  
                            E(q2)=floor(abs(Group2(j)-ClomStaticSample(i)));  
                            q2=q2+1;  
                        end  
                      end  
                  end  
                  NewCenter(2)=min(E);  
  
                  E=zeros(1,len3-1);  
                  q3=1;  
                  for j=1:len3  
                      for i=1:number  
                        if(Group3(j)~=Center(3)&&i~=j)  
                            E(q3)=floor(abs(Group3(j)-ClomStaticSample(i)));  
                            q3=q3+1;  
                        end  
                      end  
                  end  
                  NewCenter(3)=min(E);  
  
            %�ж��µ���;���ľ��������Ƿ�ͬ,��ͬ���������,����������  
            JudgeEqual=zeros(1,k);  
            for i=1:k  
                JudgeEqual=(NewCenter==Center);  
            end  
  
            S=0;  
            for i=1:k  
                if(JudgeEqual(i)==1)  
                    S=S+1;  
                end  
            end  
  
            if(S==3)  
                break;  
            end  
  
            Circulm=Circulm+1;  
     end  
     CenterSum5=zeros(time,k);           %����ÿ�γ�����kmediod�������ĵĽ��ֵ.  
     CenterSum5(i,1)=Center(1);  
     CenterSum5(i,2)=Center(2);  
     CenterSum5(i,3)=Center(3);  
end  
  
  
%����ÿ�ξ������ĵ㵽�������е�ľ���͵���Сֵ��Ϊ���ž�������  
Sum=zeros(1,time);  
for i=1:time  
    for j=1:k  
        for r=1:number-1  
            if( CenterSum5(i,j)~=ClomStaticSample(r))  
            Sum(i)=Sum(i)+CenterSum5(i,j)-ClomStaticSample(r);  
            end  
        end  
    end  
end  
  
[SumOrder, CenterEnd]=sort(Sum);%���ž������ļ�ΪCenter(CenterEnd);  
  
  
%�Դ����ݽ������յľ��ࣨ����ѡ����������ž������ģ�  
        q1=1;  
        q2=1;  
        q3=1;  
        for i=1:length(ClomStatic)  
            for j=1:3  
                EndTempDistance(i,j)=abs(ClomStatic(i)-CenterSum5(CenterEnd,j));  
            end  
            [RowMin, RowIndex]=min(EndTempDistance(i,:));  
            if(RowIndex(1)==1)  
                EndGroup1(q1)=ClomStatic(i);  
                q1=q1+1;  
            elseif(RowIndex(1)==2)  
                EndGroup2(q2)=ClomStatic(i);  
                q2=q2+1;  
            elseif(RowIndex(1)==3)  
                EndGroup3(q3)=ClomStatic(i);  
                q3=q3+1;  
            end  
        end  