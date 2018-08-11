clear;
a = [0,1,inf,inf,1;1,0,1,inf,inf;inf,1,0,1,inf;inf,inf,1,0,1;1,inf,inf,1,0];
n = length(a);
for k=1:(n-1)
    b=[1:(k-1),(k+1):n];
    kk=length(b);
    a_id=k;
    b1=[(k+1):n];
    kk1=length(b1);
    while kk>0
        for j=1:kk1
            te=a(k,a_id)+a(a_id,b1(j));
            if te<a(k,b1(j))
                a(k,b1(j))=te;
            end
        end
        miid=1;
        for j=2:kk
            if a(k,b(j))<a(k,b(miid))
                miid=j;
            end
        end
        a_id=b(miid);
        b=[b(1:(miid-1)),b((miid+1):kk)];
        kk=length(b);
        if a_id>k
            miid1=find(b1==a_id);
            b1=[b1(1:(miid-1)),b1((miid1+1):kk1)];
            kk1=length(b1);
        end
    end
    for j=(k+1):n
        a(j,k)=a(k,j);
    end
end
a;
