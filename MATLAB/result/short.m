function [d,index1,index2] = short(W)
M=inf;
p=zeros(length(W));
d=p;
d(:)=M;
for i=1:length(W)    
    pb(1:length(W))=0;
    pb(i)=1;
    index1=i;  
    index2=ones(1,length(W));
    temp=i;    
    d(i,i)=0;
    while sum(pb)<length(W)  %判断是否已经定位完成
    tb=find(pb==0);%查找pb中等于0的元素，即还未定位的那些节点
    d(i,tb)=min(d(i,tb),d(i,temp)+W(temp,tb)); %比较权值大小，将小的存入d
    tmpb=find(d(i,tb)==min(d(i,tb)));%选取最小值的节点进行下一轮的计算
    temp=tb(tmpb(1));%取相同值的第一个值作为下一轮的计算，如果有多个相同值
    pb(temp)=1;
    index1=[index1,temp];%将下一个进行的节点号接入
    index=index1(d(i,index1)==d(i,temp)-W(temp,index1));
    if length(index)>=2
        index=index(1);
    end
    index2(temp)=index;
    end
end
d;
index1;
index2;


