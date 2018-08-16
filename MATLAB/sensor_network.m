%输入样本点及其相应的类别
P = [-0.5  -0.5   0.3  -0.1  0.2   0.0   0.6  0.8;
      -0.5   0.5  -0.5   1.0  0.5  -0.9  0.8  -0.6];
T = [1 1 0 1 1 0 1 0 ];
%在坐标图上绘出样本点
plotpv(P,T);
%建立一个感知器网络
net=newp([-1 1;-1 1],1);
handle=plotpc(net.iw{1},net.b{1});
%利用样本点训练网络并绘出得到的分类线
E=1;
while (sse(E)),
   [net,Y,E]=adapt(net,P,T);
   handle=plotpc(net.iw{1},net.b{1},handle);
end;
%选择10个点来测试网络
testpoints=[-0.5  0.3 -0.9  0.4 -0.1  0.2 -0.6  0.8  0.1 -0.4;
            -0.3 -0.8 -0.4 -0.7  0.4 -0.6  0.1 -0.5 -0.5  0.3];
a=sim(net,testpoints);
%在坐标图上绘出网络的分类结果及分类线
figure;
plotpv(testpoints,a);
plotpc(net.iw{1},net.b{1},handle);
