clc;clear;

%Field Dimensions - x and y maximum (in meters)
xm=100;
ym=100;
n=200;


figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;
    S(i).G=0;
    plot(S(i).xd,S(i).yd,'o');
    hold on;
end

X=[XR;YR];
x=X';

k=4;

L = [];
L1 = 0;
while length(unique(L)) ~= k
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1))); %size(X,2)是数据集合X的数据点的数目，C是中心点的集合
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L); %-1 
        D = cumsum(sqrt(dot(D,D,1))); %将每个数据点与中心点的距离，依次累加
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1)); %find的第二个参数表示返回的索引的数目
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).')); %碉堡了，这句，将每个数据点进行分类。
    end
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
end

for i=1:length(L)
    if L(i)==1
        plot(x(i,1),x(i,2),'r*');
        plot(C(1,1),C(2,1),'r+');
        hold on;
    elseif L(i)==2
        plot(x(i,1),x(i,2),'b*');
        plot(C(1,2),C(2,2),'b+');
        hold on;
    elseif L(i)==3
        plot(x(i,1),x(i,2),'g*');
        plot(C(1,3),C(2,3),'g+');
        hold on;
    else
        plot(x(i,1),x(i,2),'k*');
        plot(C(1,4),C(2,4),'k+');
        hold on;
    end
end