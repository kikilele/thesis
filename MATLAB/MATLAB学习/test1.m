% 如果两个顶点相邻的充要条件为：它们之间的距离不大于r（r为节点的通信半径）。试选取合适的r，使得节点的平均度为12
clear;

for r=0.11:0.001:0.14
    n=200;
    p=rand(n,1)+1i*rand(n,1);
    
    %A中的第(m,n)个元素表示p中第m到第n个点的距离向量（A为NxN矩阵）,即A中每个元素代表一条线（p中m点到n点的线）
    %repmat函数是一个处理大矩阵且内容有重复时使用，其功能是以A的内容堆叠在（MxN）的矩阵B中
    %B矩阵的大小由MxN及A矩阵的内容决定，如果A是一个3x4x5的矩阵，有B = repmat(A,2,3)则最后的矩阵是6x12x5
    A=repmat(p,1,n)-repmat(p.',n,1);
    
    D=abs(A); %转换为距离，即第m到第n个点的距离
    
    %寻找所有距离小于r的连线，用在D中的位置表示，IS为行，JS为列（都是Nx1的列向量）
    [IS,JS]=find(D<r);
    plot([p(IS) p(JS)]','O-');
    print(gcf,'-dpng',num2str(r));
    %print(gcf,'-dpng',sprintf('r=%.4f.png',r));%将图形输出到文件
end