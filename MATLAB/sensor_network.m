%���������㼰����Ӧ�����
P = [-0.5  -0.5   0.3  -0.1  0.2   0.0   0.6  0.8;
      -0.5   0.5  -0.5   1.0  0.5  -0.9  0.8  -0.6];
T = [1 1 0 1 1 0 1 0 ];
%������ͼ�ϻ��������
plotpv(P,T);
%����һ����֪������
net=newp([-1 1;-1 1],1);
handle=plotpc(net.iw{1},net.b{1});
%����������ѵ�����粢����õ��ķ�����
E=1;
while (sse(E)),
   [net,Y,E]=adapt(net,P,T);
   handle=plotpc(net.iw{1},net.b{1},handle);
end;
%ѡ��10��������������
testpoints=[-0.5  0.3 -0.9  0.4 -0.1  0.2 -0.6  0.8  0.1 -0.4;
            -0.3 -0.8 -0.4 -0.7  0.4 -0.6  0.1 -0.5 -0.5  0.3];
a=sim(net,testpoints);
%������ͼ�ϻ������ķ�������������
figure;
plotpv(testpoints,a);
plotpc(net.iw{1},net.b{1},handle);
