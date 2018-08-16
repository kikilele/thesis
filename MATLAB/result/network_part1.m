function [C1,C2,b1,b2]=network_part(adj)  
C1=zeros(length(adj));
C2=zeros(length(adj));
fv=fiedlerVector(adj);
A1=find(fv>=0);
A2=find(fv<0);
a=1:length(adj);
b1 = ismember(a,A1);
b2 = ismember(a,A2);
for i=1:length(b1)
    if b1(i)~=0
        for j=1:length(b1)
            if b1(j)~=0
                C1(i,j)=adj(i,j);
            end
        end
    end
end

for i=1:length(b2)
    if b2(i)~=0
        for j=1:length(b2)
            if b2(j)~=0
                C2(i,j)=adj(i,j);
            end
        end
    end
end

C1(all(C1==0,2),:) =[];
C1(:,all(C1==0,1)) =[];
C2(all(C2==0,2),:) =[];
C2(:,all(C2==0,1)) =[];

C1;
C2;
b1;
b2;