% �б����ʱ��ͨ���漰��������������֮��ľ��룬��Ԫͳ��ѧ�������ж��־�����㹫ʽ��
% MATLAB�����ж�Ӧ�������ɷ���ֱ�ӵ��ü��㡣
% ���뺯���У�pdist, pdist2, mahal, squareform, mdscale, cmdscale
% ��Ҫ����pdist2 ,�����ɲο�matlab help
%  
% D = pdist2(X,Y)
% D = pdist2(X,Y,distance)
% D = pdist2(X,Y,'minkowski',P)
% D = pdist2(X,Y,'mahalanobis',C)
% D = pdist2(X,Y,distance,'Smallest',K)
% D = pdist2(X,Y,distance,'Largest',K)
% [D,I] = pdist2(X,Y,distance,'Smallest',K)
% [D,I] = pdist2(X,Y,distance,'Largest',K)

% 2�ּ��㷽ʽ��һ��ֱ������pdist���㣬��һ�ְ���ʽֱ�Ӽ���,����һ��Ա������
clc;clear;
x = rand(4,3);
y = rand(1,3);
for i =1:size(x,1)
    for j =1:size(y,1)
        a = x(i,:); b=y(j,:);
       
%         Euclidean distance
        d1(i,j)=sqrt((a-b)*(a-b)');
       
%         Standardized Euclidean distance
        V = diag(1./std(x).^2);
        d2(i,j)=sqrt((a-b)*V*(a-b)');
       
%         Mahalanobis distance
        C = cov(x);
        d3(i,j)=sqrt((a-b)*pinv(C)*(a-b)');
       
%         City block metric
        d4(i,j)=sum(abs(a-b));
       
%         Minkowski metric
        p=3;
        d5(i,j)=(sum(abs(a-b).^p))^(1/p);
       
%         Chebychev distance
        d6(i,j)=max(abs(a-b));
       
%         Cosine distance
        d7(i,j)=1-(a*b')/sqrt(a*a'*b*b');
       
%         Correlation distance
        ac = a-mean(a); bc = b-mean(b);       
        d8(i,j)=1- ac*bc'/(sqrt(sum(ac.^2))*sqrt(sum(bc.^2)));
    end
end

md1 = pdist2(x,y,'Euclidean');
md2 = pdist2(x,y,'seuclidean');
md3 = pdist2(x,y,'mahalanobis');
md4 = pdist2(x,y,'cityblock');
md5 = pdist2(x,y,'minkowski',p);
md6 = pdist2(x,y,'chebychev');
md7 = pdist2(x,y,'cosine');
md8 = pdist2(x,y,'correlation');
md9 = pdist2(x,y,'hamming');
md10 = pdist2(x,y,'jaccard');
md11 = pdist2(x,y,'spearman');
D1=[d1,md1],D2=[d2,md2],D3=[d3,md3]
D4=[d4,md4],D5=[d5,md5],D6=[d6,md6]
D7=[d7,md7],D8=[d8,md8]
