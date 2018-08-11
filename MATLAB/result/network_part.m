function [b1,b2]=network_part(adj)  
fv=fiedlerVector(adj);
A1=find(fv>=0);
A2=find(fv<0);
b1=[];
b2=[];
for i=1:length(adj)
   if ismember(i,A1)
       b1=[b1,i];
   end
end
for i=1:length(adj)
   if ismember(i,A2)
       b2=[b2,i];
   end
end
b1;
b2;
