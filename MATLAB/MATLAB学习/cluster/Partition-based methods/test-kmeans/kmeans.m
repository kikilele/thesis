clc;  
clear;  
  
ClomStatic=[1,2,3,25,26,27,53,54,55];  
len=length(ClomStatic);%求向量ClomStatic的长度  
  
k=3; %给定的类别数目  
  
%产生三个随机整数，随机聚类中心  
p=randperm(len);  
Temp=p(1:k);  
Center=zeros(1,k);  
for i=1:k  
    Center(i)=ClomStatic(Temp(i));  
end  
  
  
  
%计算除聚类中心外的样本数据到聚类中心的距离,然后进行聚类  
TempDistance=zeros(len,3);  
 while 1  
      
    Circulm=1;  
      
    p1=1;  
    p2=1;  
    p3=1;  
      
    JudgeEqual=zeros(1,k);  
    if(Circulm~=1)  
        clear Group1 Group2 Group3;     
    end  
    for i=1:len  
        for j=1:3  
            TempDistance(i,j)=abs(ClomStatic(i)-Center(j));  
        end  
        [RowMin, RowIndex]=min(TempDistance(i,:));  
        if(RowIndex==1)  
            Group1(p1)=ClomStatic(i);  
            p1=p1+1;  
        elseif(RowIndex==2)  
            Group2(p2)=ClomStatic(i);  
            p2=p2+1;  
        elseif(RowIndex==3)  
            Group3(p3)=ClomStatic(i);  
            p3=p3+1;  
        end  
    end  
       
        len1=length(Group1);  
        len2=length(Group2);  
        len3=length(Group3);  
          
          
        %计算Group1,Group2,Group3的均值  
        MeanGroup1=mean(Group1);  
        MeanGroup2=mean(Group2);  
        MeanGroup3=mean(Group3);  
  
          
        %分别计算每个类中距离均值最近的一个点作为新的聚类中心  
        AbsGroup1=zeros(1,len1);  
        for t=1:len1  
            AbsGroup1(t)=floor(abs(Group1(t)-MeanGroup1));  
        end  
        [MaxAbsGroup1, MaxAbsGroup1Index]=min(AbsGroup1);  
        NewCenter(1)=Group1(MaxAbsGroup1Index);  
        clear AbsGroup1;  
  
        AbsGroup2=zeros(1,len2);  
        for t=1:len2  
            AbsGroup2(t)=floor(abs(Group2(t)-MeanGroup2));  
        end  
        [MaxAbsGroup2, MaxAbsGroup2Index]=min(AbsGroup2);  
        NewCenter(2)=Group2(MaxAbsGroup2Index);  
        clear AbsGroup2;  
            
        AbsGroup3=zeros(1,len3);  
        for t=1:len3  
            AbsGroup3(t)=floor(abs(Group3(t)-MeanGroup3));  
        end  
        [MaxAbsGroup3, MaxAbsGroup3Index]=min(AbsGroup3);  
        NewCenter(3)=Group3(MaxAbsGroup2Index);  
        clear AbsGroup3;  
          
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