clear
clc
%row\col：图片大小
%img：图
%dm：扩展边长
%blocksize：分割方式
%rowblocks\colblocks：图片分割的块数
%SAD：绝对差值和
%mvx\mvy：运动矢量
%diff：数据流动
%pre：预测下一帧图像

%显示图片
img1=imread('1.jpg'); 
img2=imread('2.jpg');
%figure(1)
%imshow(img1)
%title('第一帧原图')
%figure(2)
%imshow(img2)
%title('第二帧原图')
img1=rgb2gray(img1);
img2=rgb2gray(img2);
I1=double(img1);
I2=double(img2);
%figure(3)
%imshow(img1)
%title('第一帧灰度图')
%figure(4)
%imshow(img2)
%title('第二帧灰度图')
%I3=I2-I1;
%figure(5)
%imshow(I3)
%title('帧间差值')
[row col]=size(I1);
dm=7;
II=zeros(row+2*dm,col+2*dm);
II(dm+1:dm+row,dm+1:dm+col)=I1;
for i=1:dm
    II(i,dm+1:dm+col)=II(dm+1,dm+1:dm+col);
    II(row+dm+i,dm+1:dm+col)=II(dm+row,dm+1:dm+col);
end
for j=1:dm
    II(1:row+2*dm,j)=II(1:row+2*dm,dm+1);
    II(1:row+2*dm,col+dm+j)=II(1:row+2*dm,dm+col);
end

%全搜索算法
 blocksize=16;
 rowblocks =ceil(row/blocksize);
 colblocks =ceil(col/blocksize);
 A=100000000000000000000;         
 %mvx=ones(rowblocks,colblocks);           
 %mvy=ones(rowblocks,colblocks); 
 mvx=zeros(rowblocks,colblocks);           
 mvy=zeros(rowblocks,colblocks); 
 diff=zeros(row,col);         
 tic
 for x=0:(rowblocks-1)         
     rownum=x*blocksize;
     for y=0:(colblocks-1)       
         colnum=y*blocksize;
         for p=-dm:dm
             for q=-dm:dm      
                 SAD=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+p+1:rownum+dm+p+blocksize,colnum+dm+q+1:colnum+dm+q+blocksize)))); 
                 if SAD<A
                    A=SAD;
                    mvx(x+1,y+1)=p;
                    mvy(x+1,y+1)=q;
                 end   
             end
         end
         A=100000000000000000000;
         for mx=1:blocksize
             for ny=1:blocksize
                 diff(rownum+mx,colnum+ny)=I2(rownum+mx,colnum+ny)-II(rownum+mx+dm+mvx(x+1,y+1),colnum+ny+dm+mvy(x+1,y+1));
             end
         end
         %pre(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)=II(rownum+dm+mvx(x+1,y+1)+1:rownum+dm+mvx(x+1,y+1)+blocksize,colnum+dm+mvy(x+1,y+1)+1:colnum+dm+mvy(x+1,y+1)+blocksize)+diff(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize);
         pre1(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)=II(rownum+dm+mvx(x+1,y+1)+1:rownum+dm+mvx(x+1,y+1)+blocksize,colnum+dm+mvy(x+1,y+1)+1:colnum+dm+mvy(x+1,y+1)+blocksize);
     end
 end
toc
%figure(6)
%imshow(diff,[]);
%title('匹配后的帧间差值')
%figure(7)
%imshow(pre1,[]);
%title('匹配后的第二帧图像');
%figure(8)
%imshow(pre,[]);
%title('恢复后的第二帧图像');

%绘制运动矢量图
%figure(9)
 %for i=1:16
     %for j=1:16
         %quiver(i,j,mvx(i,j)/16,mvy(i,j)/16);hold on
     %end
 %end
%quiver(1:16,1:16,mvy,mvx)
%quiver(1:colblocks,1:rowblocks,mvy,mvx)
%title('运动矢量图')
%grid on
                   
%PSNR
%MSE=sum(sum((I2-pre1).*(I2-pre1)))/256; 
%MSE=sum(sum((I2-pre1).*(I2-pre1)))/(rowblocks*colblocks); 
err=0;
for i=1:row
    for j=1:col
        err=err+(I2(i,j)-pre1(i,j))*(I2(i,j)-pre1(i,j));
    end
end
MSE=err/(row*col);
PSNR=10*log10(255*255/MSE);
fprintf('the PSNR of predicted image:%6.3f dB\n',PSNR)