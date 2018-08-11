clear;
xm=100;
ym=100;
n=100;
C=cell(1,20);%´æfv
B=cell(1,20);
k=2;
b1=[];
b2=[];
cluster=2;

figure(1);
for i=1:1:n
    S(i).xd=rand(1,1)*xm;%×ø±ê
    XR(i)=S(i).xd;
    S(i).yd=rand(1,1)*ym;
    YR(i)=S(i).yd;    
    ID(i)=i;
    plot(S(i).xd,S(i).yd,'ko');
    hold on; 
end

for i=1:n
    for j=1:n
        D(i,j)=sqrt((XR(i)-XR(j))^2 + (YR(i)-YR(j))^2 );
        if D(i,j)>50
            D(i,j)=0;   
        elseif i==j 
            D(i,j)=0;
        else D(i,j)=1;
        end
    end
end

%[b1,b2]=network_part(D);
fv=fiedlerVector(D);
A1=find(fv>=0);
A2=find(fv<0);
for i=1:length(D)
   if ismember(i,A1)
       b1=[b1,i];
   end
end
for i=1:length(D)
   if ismember(i,A2)
       b2=[b2,i];
   end
end
B{1}=b1;
B{2}=b2;
for i=1:length(B{1})
    for j=1:length(B{1})
        C1(i,j)=sqrt((XR(B{1}(1,i))-XR(B{1}(1,j)))^2 + (YR(B{1}(1,i))-YR(B{1}(1,j)))^2 );
        if C1(i,j)>50
            C1(i,j)=0;   
        elseif i==j 
            C1(i,j)=0;
        else C1(i,j)=1;
        end
    end
end
C{1}=C1;
s1=graphSpectrum(C1);
for i=1:length(B{2})
    for j=1:length(B{2})
        C2(i,j)=sqrt((S(B{2}(1,i)).xd-S(B{2}(1,j)).xd)^2 + (S(B{2}(1,i)).yd-S(B{2}(1,j)).yd)^2 );
        if C2(i,j)>50
            C2(i,j)=0;   
        elseif i==j 
            C2(i,j)=0;
        else C2(i,j)=1;
        end
    end
end
C{2}=C2;
s2=graphSpectrum(C2);
P=[s1,s2];

figure(2);
for i=1:length(B{1})
    plot(XR(B{1}(1,i)),YR(B{1}(1,i)),'ko');
    hold on;
end
for i=1:length(B{2})
    plot(XR(B{2}(1,i)),YR(B{2}(1,i)),'ro');
    hold on;
end


while (cluster<4)    
    cluster=cluster+1;
    m=length(P);
    P_min=min(P);
    for e=1:m  
        if P(e)==P_min
            P(e)=inf;
            b1=[];
            b2=[];
            A1=[];
            A2=[];
            %[b1,b2]=network_part(C{e});
            fv=fiedlerVector(C{e});
            A1=find(fv>=0);
            A2=find(fv<0);
            for i=1:length(C{e})
                if ismember(i,A1)
                    b1=[b1,B{e}(1,i)];
                end
            end
            for i=1:length(C{e})
                if ismember(i,A2)
                    b2=[b2,B{e}(1,i)];
                end
            end
            C{e}=[];
            B{e}=[];
            k=k+1;
            B{k}=b1;
            C1=[];
            for i=1:length(b1)
                for j=1:length(b1)                    
                    C1(i,j)=sqrt((S(b1(i)).xd-S(b1(j)).xd)^2 + (S(b1(i)).yd-S(b1(j)).yd)^2 );
                    if C1(i,j)>50
                        C1(i,j)=0;
                    elseif i==j
                        C1(i,j)=0;
                    else C1(i,j)=1;
                    end
                end
            end
            C{k}=C1;
            s1=graphSpectrum(C1);
            k=k+1;            
            B{k}=b2;
            C2=[];
            for i=1:length(b2)
                for j=1:length(b2)                    
                    C2(i,j)=sqrt((S(b2(i)).xd-S(b2(j)).xd)^2 + (S(b2(i)).yd-S(b2(j)).yd)^2 );
                    if C2(i,j)>50
                        C2(i,j)=0;
                    elseif i==j
                        C2(i,j)=0;
                    else C2(i,j)=1;
                    end
                end
            end
            C{k}=C2;
            s2=graphSpectrum(C2);
            P=[P,s1,s2];            
        end       
    end
    figure(cluster);
    for i=1:length(B)
        if ~isempty(B{i})
            switch i
                case {1}
                    for j=1:length(B{1})
                        plot(XR(B{1}(1,j)),YR(B{1}(1,j)),'ko');
                        hold on;
                    end
                case {2}
                    for j=1:length(B{2})
                        plot(XR(B{2}(1,j)),YR(B{2}(1,j)),'ro');
                        hold on;
                    end
                case {3}
                    for j=1:length(B{3})
                        plot(XR(B{3}(1,j)),YR(B{3}(1,j)),'go');
                        hold on;
                    end    
                 case {4}
                    for j=1:length(B{4})
                        plot(XR(B{4}(1,j)),YR(B{4}(1,j)),'bo');
                        hold on;
                    end 
                 case {5}
                    for j=1:length(B{5})
                        plot(XR(B{5}(1,j)),YR(B{5}(1,j)),'co');
                        hold on;
                    end 
                 case {6}
                    for j=1:length(B{6})
                        plot(XR(B{6}(1,j)),YR(B{6}(1,j)),'yo');
                        hold on;
                    end    
            end
        end
    end
end




for i=1:length(C)
    if ~isempty(C{i})        
        PageRank(C{i},0.001)
    end
end
    
