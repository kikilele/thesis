%两点间最短路的Dijksrta算法
function [d, index1, index2] = Dijkf(a)
% a 表示图的权值矩阵
% d 表示所求最短路的权和
% index1表示标号顶点顺序
% index2表示标号顶点索引
% 只能计算从单一原点开始

%参数初始化
M=inf;
pb(1:length(a))=0;%标号的信息，为1表示该节点已经标号，初始化全为0
pb(1)=1;
index1 =1;%表示从哪一个开始
index2=ones(1,length(a));%表示到此节点的上一个节点号,初始化为index1的值
d(1:length(a))=M;%初始化为无穷
d(1)=0;
temp=1;
while sum(pb)<length(a)  %判断是否已经定位完成
    tb=find(pb==0);%查找pb中等于0的元素，即还未定位的那些节点
    d(tb)=min(d(tb),d(temp)+a(temp,tb)); %比较权值大小，将小的存入d
    tmpb=find(d(tb)==min(d(tb)));%选取最小值的节点进行下一轮的计算
    temp=tb(tmpb(1));%取相同值的第一个值作为下一轮的计算，如果有多个相同值
    pb(temp)=1;
    index1=[index1,temp];%将下一个进行的节点号接入
    index=index1(d(index1)==d(temp)-a(temp,index1));
    if length(index)>=2
        index=index(1);
    end
    index2(temp)=index;
end
d;
index1;
index2;











