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
    while sum(pb)<length(W)  %�ж��Ƿ��Ѿ���λ���
    tb=find(pb==0);%����pb�е���0��Ԫ�أ�����δ��λ����Щ�ڵ�
    d(i,tb)=min(d(i,tb),d(i,temp)+W(temp,tb)); %�Ƚ�Ȩֵ��С����С�Ĵ���d
    tmpb=find(d(i,tb)==min(d(i,tb)));%ѡȡ��Сֵ�Ľڵ������һ�ֵļ���
    temp=tb(tmpb(1));%ȡ��ֵͬ�ĵ�һ��ֵ��Ϊ��һ�ֵļ��㣬����ж����ֵͬ
    pb(temp)=1;
    index1=[index1,temp];%����һ�����еĽڵ�Ž���
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


