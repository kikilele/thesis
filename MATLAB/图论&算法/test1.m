clear;
a = [0,1,inf,inf,1;1,0,1,inf,inf;inf,1,0,1,inf;inf,inf,1,0,1;1,inf,inf,1,0];
M=inf;
p=zeros(length(a));
%pb=p;
d=p;
%index1=p;
d(:)=M;
%index2=ones(length(a));
%pb(1:length(a))=0;
%pb(1)=1;
%index1 =1;
%index2=ones(1,length(a));
%d(1:length(a))=M;
%d(1)=0;
%temp=1;
for i=1:length(a)    
    pb(1:length(a))=0;
    pb(i)=1;
    index1=i;  
    index2=ones(1,length(a));
    temp=i;    
    d(i,i)=0;
    while sum(pb)<length(a)  %�ж��Ƿ��Ѿ���λ���
    tb=find(pb==0);%����pb�е���0��Ԫ�أ�����δ��λ����Щ�ڵ�
    d(i,tb)=min(d(i,tb),d(i,temp)+a(temp,tb)); %�Ƚ�Ȩֵ��С����С�Ĵ���d
    tmpb=find(d(i,tb)==min(d(i,tb)));%ѡȡ��Сֵ�Ľڵ������һ�ֵļ���
    temp=tb(tmpb(1));%ȡ��ֵͬ�ĵ�һ��ֵ��Ϊ��һ�ֵļ��㣬����ж����ֵͬ
    pb(temp)=1;
    index1=[index1,temp];%����һ�����еĽڵ�Ž���
    index=index1(d(i,index1)==d(i,temp)-a(temp,index1));
    if length(index)>=2
        index=index(1);
    end
    index2(temp)=index;
    end
end
d;