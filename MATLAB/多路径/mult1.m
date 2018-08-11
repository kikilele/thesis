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
img1=imread('image13.jpg'); 
img2=imread('image14.jpg');
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

%多路径搜索法
 blocksize=16;
 rowblocks=ceil(row/blocksize);
 colblocks=ceil(col/blocksize);        
 mvx=zeros(rowblocks,colblocks);           
 mvy=zeros(rowblocks,colblocks); 
 diff=zeros(row,col);
 tic  
for x=0:(rowblocks-1)       
     rownum=x*blocksize;
     for y=0:(colblocks-1)        
         colnum=y*blocksize;
         M=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+1:rownum+dm+blocksize,colnum+dm+1:colnum+dm+blocksize))));
         MIN=M;                   %第一步：得到原点SAD，定为MIN
         pointx=zeros(8,10);
         pointy=zeros(8,10);
         k=1;
         for p1=-1:1              %第二步：搜索原点周围的8个点
             for q1=-1:1   
                 if abs(p1)~=0 | abs(q1)~=0
                    SAD=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+p1+1:rownum+dm+p1+blocksize,colnum+dm+q1+1:colnum+dm+q1+blocksize))));
                    if SAD<M
                       M=SAD;
                    end
                    if ~(SAD<MIN)        %原点为最佳匹配点
                        %mvx(x+1,y+1)=0;
                        %mvy(x+1,y+1)=0;
                        break;
                    else
                        %for j=1:10
                        pointx(1,k)=p1;
                        pointy(1,k)=q1;
                        k=k+1;
                        %end
                    end
                 end
             end
             if ~(SAD<MIN)
                 break;
             end
         end
         MIN=M;
         if ~(pointx==zeros(8,10) & pointy==zeros(8,10)) %第三步：对锚点进行搜索
            for i=2:8    
               for j=1:10
                   k=1;
                   for p1=pointx(i-1,j)-1:pointx(i-1,j)+1
                       for q1=pointy(i-1,j)-1:pointy(i-1,j)+1
                           if abs(p1)>=i-1 | abs(q1)>=i-1
                              SAD=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+p1+1:rownum+dm+p1+blocksize,colnum+dm+q1+1:colnum+dm+q1+blocksize))));  
                              if SAD<M
                                 M=SAD;
                              end
                              if ~(SAD<MIN)
                                 break; 
                              else
                                 %for k=1:10
                                  p1=pointx(i,k);
                                  q1=pointy(i,k);
                                  k=k+1;
                                 %end
                              end
                           end
                       end
                       if ~(SAD<MIN)
                       break;
                       end
                   end
               end
            end
         end
         for i=1:8
             for j=1:10
                 if ~(pointx(i,j)==0 & pointy(i,j)==0)
                    mvx(x+1,y+1)=pointx(i,j);
                    mvy(x+1,y+1)=pointy(i,j);
                 end
             end
         end
         for mx=1:blocksize
             for ny=1:blocksize
                 diff(rownum+mx,colnum+ny)=I2(rownum+mx,colnum+ny)-II(rownum+mx+dm+mvx(x+1,y+1),colnum+ny+dm+mvy(x+1,y+1));
             end
         end
         pre(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)=II(rownum+dm+mvx(x+1,y+1)+1:rownum+dm+mvx(x+1,y+1)+blocksize,colnum+dm+mvy(x+1,y+1)+1:colnum+dm+mvy(x+1,y+1)+blocksize);
     end
end
toc
%figure(6)
%imshow(diff,[]);
%title('匹配后的帧间差值')
%figure(7)
%imshow(pre,[]);
%title('匹配后的第二帧图像');

%绘制运动矢量图 
%figure(8)
%quiver(1:colblocks,1:rowblocks,mvy,mvx)
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