clear;
W=[0,1,1,3;1,0,2,inf;1,2,0,2;3,inf,2,0];
n=length(W);
U=W;
m=1;
while m<=n
    for i=1:n
        for j=1:n
            if U(i,j)>U(i,m)+U(m,j)
                U(i,j)=U(i,m)+U(m,j);
            end
        end
    end
    m=m+1;
end
u=U(1,n);
P1=zeros(1,n);
k=1;
P1(k)=n;
V=ones(1,n)*inf;
kk=n;
while kk~=1
    for i=1:n
        V(1,i)=U(1,kk)-W(i,kk);
        if V(1,i)==U(1,i)
            P1(k+1)=i;
            kk=i;
            k=k+1;
        end
    end
end
k=1;
wrow=find(P1~=0);
for j=length(wrow):(-1):1
    P(k)=P1(wrow(j));
    k=k+1;
end
P;%默认是输出从1到n的最短路径
u;%是最短路径的权
