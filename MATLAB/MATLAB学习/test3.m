clear;

n=200;
p=rand(n,1)+1i*rand(n,1);


% for i=1:length(p)
%     text(real(p(i)),imag(p(i)),num2str(i),'Fontsize',14,'color','r');    
% end

A=repmat(p,1,n)-repmat(p.',n,1);
D=abs(A);
save figuredate.mat p A D;

for r=0.2:0.001:0.3
    [IS,JS]=find(D<r);
    sizeIS=size(IS);
    avg=(sizeIS(1)-n)/2/n;
    fprintf('r=%.4f,Æ½¾ù¶È=%.4f\n',r,avg);
    if avg>=12
        plot([p(IS) p(JS)]','O-');
        break;
    end
end