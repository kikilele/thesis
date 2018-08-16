clc;  
clear;  
  
load Data3.mat;  
  
k=3; %给定的类别数目  
  
time=5;%time为抽样的次数  
number=30;%number为抽样本个数  
for T=1:time  
    ClomStaticSample=zeros(1,number);  
    ClomStaticSample=randsample(ClomStatic,number);   %ClomStaticSample就是样本数据  
                                                      %接下来对样本数据使用kmediod算法进行聚类  
                                                        
    %产生三个随机整数，随机聚类中心  
    p=randperm(number);  
    Temp=p(1:k);  
    Center=zeros(1,k);  
    for j=1:k  
        Center(j)=ClomStaticSample(Temp(j));  
    end  
    [ClomStaticSample]=sort(ClomStaticSample);  
      
    TempDistance=zeros(number,3);           %暂存差值  
      
     while 1  
        Circulm=1;                          %循环控制  
  
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
              
  
                  %分别计算每个类中除开类中心的点到其他所有点的距离和E,E最小时为该类新的聚类中心.  
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
  
            %判断新的类和旧类的聚类中心是否不同,不同则继续聚类,否则聚类结束  
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
     CenterSum5=zeros(time,k);           %保存每次抽样后kmediod聚类中心的结果值.  
     CenterSum5(i,1)=Center(1);  
     CenterSum5(i,2)=Center(2);  
     CenterSum5(i,3)=Center(3);  
end  
  
  
%计算每次聚类中心点到其他所有点的距离和的最小值即为最优聚类中心  
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
  
[SumOrder, CenterEnd]=sort(Sum);%最优聚类中心即为Center(CenterEnd);  
  
  
%对大数据进行最终的聚类（按照选择出来的最优聚类中心）  
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