%��������·��Dijksrta�㷨
function [d, index1, index2] = Dijkf(a)
% a ��ʾͼ��Ȩֵ����
% d ��ʾ�������·��Ȩ��
% index1��ʾ��Ŷ���˳��
% index2��ʾ��Ŷ�������
% ֻ�ܼ���ӵ�һԭ�㿪ʼ

%������ʼ��
M=inf;
pb(1:length(a))=0;%��ŵ���Ϣ��Ϊ1��ʾ�ýڵ��Ѿ���ţ���ʼ��ȫΪ0
pb(1)=1;
index1 =1;%��ʾ����һ����ʼ
index2=ones(1,length(a));%��ʾ���˽ڵ����һ���ڵ��,��ʼ��Ϊindex1��ֵ
d(1:length(a))=M;%��ʼ��Ϊ����
d(1)=0;
temp=1;
while sum(pb)<length(a)  %�ж��Ƿ��Ѿ���λ���
    tb=find(pb==0);%����pb�е���0��Ԫ�أ�����δ��λ����Щ�ڵ�
    d(tb)=min(d(tb),d(temp)+a(temp,tb)); %�Ƚ�Ȩֵ��С����С�Ĵ���d
    tmpb=find(d(tb)==min(d(tb)));%ѡȡ��Сֵ�Ľڵ������һ�ֵļ���
    temp=tb(tmpb(1));%ȡ��ֵͬ�ĵ�һ��ֵ��Ϊ��һ�ֵļ��㣬����ж����ֵͬ
    pb(temp)=1;
    index1=[index1,temp];%����һ�����еĽڵ�Ž���
    index=index1(d(index1)==d(temp)-a(temp,index1));
    if length(index)>=2
        index=index(1);
    end
    index2(temp)=index;
end
d;
index1;
index2;











