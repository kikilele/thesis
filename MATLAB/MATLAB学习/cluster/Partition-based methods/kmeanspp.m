
function [L,C] = kmeanspp(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeanspp(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.
%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
L = [];
L1 = 0;
while length(unique(L)) ~= k
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1))); %size(X,2)�����ݼ���X�����ݵ����Ŀ��C�����ĵ�ļ���
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L); %-1 
        D = cumsum(sqrt(dot(D,D,1))); %��ÿ�����ݵ������ĵ�ľ��룬�����ۼ�
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1)); %find�ĵڶ���������ʾ���ص���������Ŀ
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).')); %�ﱤ�ˣ���䣬��ÿ�����ݵ���з��ࡣ
    end
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
end