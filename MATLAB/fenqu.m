clear;
A = [0 1 0 0 1 0;1 0 1 0 1 0;0 1 0 1 0 0;0 0 1 0 1 1;1 1 0 1 0 0;0 0 0 1 0 0];
k=1;
while (k<=10)
    L =diag(sum(A))-A;
    [V,D]=eig(L);
    X2 = V(:,2);
    A1=find(X2>=0);
    A2=find(X2<0);
    a=1:length(A);
    
    
    
    
    k=k+1;
end





b1 = ismember(a,A1);
b2 = ismember(a,A2);
C1=zeros(6);
C2=zeros(6);
for i=1:length(b1)
    if b1(i)~=0
        for j=1:length(b1)
            if b1(j)~=0
                C1(i,j)=A(i,j);
            end
        end
    end
end
C1(all(C1==0,2),:) =[];
C1(:,all(C1==0,1)) =[];
for i=1:length(b2)
    if b2(i)~=0
        for j=1:length(b2)
            if b2(j)~=0
                C2(i,j)=A(i,j);
            end
        end
    end
end
C2(all(C2==0,2),:) =[];
C2(:,all(C2==0,1)) =[];
