clear
clc

N=1000;
P(1:6)=log([2.3 2.5 2.8 2.89 3.2 3.3]);

for t=6:N-1
    temp=0;
    for i=1:5
        temp=temp+abs(P(t-i+1)-P(t-i));
    end
    judge1=(1/5)*temp;
    judge2=abs(P(t)-P(t-1));
    if (judge1<0.0225) &&(judge2<0.05) 
        Xt=0.025;
    else
        Xt=0.05;
    end
    Dt=Xt-(2/Xt)*(P(t)-P(t-1))^2;
    P(t+1)=P(t)+Dt;
end
P=P';
figure
plot(P)



