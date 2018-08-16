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
img1=imread('33.jpg'); 
img2=imread('34.jpg');
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

%����������
 blocksize=16;
 rowblocks=ceil(row/blocksize);
 colblocks=ceil(col/blocksize);
 A=99999999999999999999;         
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
         for p1=-4:4:4              %��һ��
             for q1=-4:4:4 
                 SAD=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+p1+1:rownum+dm+p1+blocksize,colnum+dm+q1+1:colnum+dm+q1+blocksize)))); 
                 if SAD<A
                    A=SAD;
                    mvx(x+1,y+1)=p1;
                    mvy(x+1,y+1)=q1;
                 end   
             end
         end
         p1=mvx(x+1,y+1);
         q1=mvy(x+1,y+1);
         for p2=p1-2:2:p1+2         %�ڶ���
             for q2=q1-2:2:q1+2
                 if p2~=p1 | q2~=q1
                    SAD=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+p2+1:rownum+dm+p2+blocksize,colnum+dm+q2+1:colnum+dm+q2+blocksize))));
                    if SAD<A
                       A=SAD;
                       mvx(x+1,y+1)=p2;
                       mvy(x+1,y+1)=q2;
                    end 
                 end
             end
         end
         p2=mvx(x+1,y+1);
         q2=mvy(x+1,y+1);
         for p3=p2-1:1:p2+1        %������
             for q3=q2-1:1:q2+1
                 if p3~=p2 | q3~=q2 
                   SAD=sum(sum(abs(I2(rownum+1:rownum+blocksize,colnum+1:colnum+blocksize)-II(rownum+dm+p3+1:rownum+dm+p3+blocksize,colnum+dm+q3+1:colnum+dm+q3+blocksize))));
                   if SAD<A
                      A=SAD;
                      mvx(x+1,y+1)=p3;
                      mvy(x+1,y+1)=q3;
                   end 
                 end
             end
         end
         A=999999999999999999;
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
         %quiver(i,j,mvx(i,j)/16,mvy(i,j)/16); hold on;
     %end
 %end
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