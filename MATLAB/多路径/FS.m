clear
clc
%row\col��ͼƬ��С
%img��ͼ
%dm����չ�߳�
%blocksize���ָʽ
%rowblocks\colblocks��ͼƬ�ָ�Ŀ���
%SAD�����Բ�ֵ��
%mvx\mvy���˶�ʸ��
%diff����������
%pre��Ԥ����һ֡ͼ��

%��ʾͼƬ
img1=imread('1.jpg'); 
img2=imread('2.jpg');
%figure(1)
%imshow(img1)
%title('��һ֡ԭͼ')
%figure(2)
%imshow(img2)
%title('�ڶ�֡ԭͼ')
img1=rgb2gray(img1);
img2=rgb2gray(img2);
I1=double(img1);
I2=double(img2);
%figure(3)
%imshow(img1)
%title('��һ֡�Ҷ�ͼ')
%figure(4)
%imshow(img2)
%title('�ڶ�֡�Ҷ�ͼ')
%I3=I2-I1;
%figure(5)
%imshow(I3)
%title('֡���ֵ')
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

%ȫ�����㷨
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
%title('ƥ����֡���ֵ')
%figure(7)
%imshow(pre1,[]);
%title('ƥ���ĵڶ�֡ͼ��');
%figure(8)
%imshow(pre,[]);
%title('�ָ���ĵڶ�֡ͼ��');

%�����˶�ʸ��ͼ
%figure(9)
 %for i=1:16
     %for j=1:16
         %quiver(i,j,mvx(i,j)/16,mvy(i,j)/16);hold on
     %end
 %end
%quiver(1:16,1:16,mvy,mvx)
%quiver(1:colblocks,1:rowblocks,mvy,mvx)
%title('�˶�ʸ��ͼ')
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