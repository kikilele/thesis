Ant=100;
ECHO=50;
start1=-1;end1=1;start2=-1;end2=1;tcl=0.05;
f='cos(2*pi.*x).*cos(2*pi.*y).*exp(-((x.^2+y.^2)/10))';
[x,y]=meshgrid(start1:tcl:end1,start2:tcl:end2,start3:tcl:end3,start4:tcl:end4);vxp=x;vyp=y;

vzp=eval(f);
for i=1:Ant
    X(i,1)=(start1+(end1-start1)*rand(1));
    X(i,2)=(start2+(end2-start2)*rand(1));
    
    T0(i)=cos(2*pi.*X(i,1)).*cos(2*pi.*X(i,2)).*exp(-((X(i,1).^2+X(i,2).^2)/10));
end;

figure(1);mesh(vxp,vyp,vzp);hold on;
plot3(X(:,1),X(:,2),T0,'k*')
hold on;grid on;title('(a)');xlabel('x');ylabel('y');zlabel('f(x,y)');
for Echo=1:ECHO
    P0=0.2;P=0.8;
    lamda=1/Echo;
    [T_Best(Echo),BestIndex]=max(T0);
    for j_g=1:Ant
        r=T0(BestIndex)-T0(j_g);
        Prob(Echo,j_g)=r/T0(BestIndex);
    end;
     for j_g_tr=1:Ant
         if Prob(Echo,j_g_tr)<P0
           temp1=X(j_g_tr,1)+(2*rand(1)-1)*lamda;
           temp2=X(j_g_tr,2)+(2*rand(1)-1)*lamda;
         else
             temp1=X(j_g_tr,1)+(end1-start1)*(rand(1)-0.5);
             temp2=X(j_g_tr,2)+(end2-start2)*(rand(1)-0.5);
         end
         if temp1<start1
            temp1=start1;
         end;
         if temp1>end1
             temp1=end1;
         end;
         if temp2<start2
            temp2=start2;
         end;
         if temp2>end2
             temp2=end2;
         end;
         if cos(2*pi.*temp1).*cos(2*pi.*temp2).*exp(-((temp1.^2+temp2.^2)/10))>cos(2*pi.*X(j_g_tr,1)).*cos(2*pi.*X(j_g_tr,2)).*exp(-((X(j_g_tr,1).^2+X(j_g_tr,2).^2)/10))
             X(j_g_tr,1)=temp1;X(j_g_tr,2)=temp2;
         end;
     end;
     for t_t=1:Ant
         T0(t_t)=(1-P)*T0(t_t)+cos(2*pi.*X(t_t,1)).*cos(2*pi.*X(t_t,2)).*exp(-((X(t_t,1).^2+X(t_t,2).^2)/10));
     end;
     [c_iter,i_iter]=max(T0);
     maxpiont_iter=[X(i_iter,1),X(i_iter,2)];
     maxvalue_iter=cos(2*pi.*X(i_iter,1)).*cos(2*pi.*X(i_iter,2)).*exp(-((X(i_iter,1).^2+X(i_iter,2).^2)/10));
     max_local(Echo)=maxvalue_iter;
     if Echo>=2
         if max_local(Echo)>max_global(Echo-1)
             max_global(Echo)=max_local(Echo);
         else
             max_global(Echo)=max_global(Echo-1);
         end
     else
         max_global(Echo)=maxvalue_iter;
     end;
end;%ECHO循环结束
figure(2);mesh(vxp,vyp,vzp);hold on;
x=X(:,1);y=X(:,2);
plot3(x,y,eval(f),'k*')
hold on;grid on;title('蚂蚁的最终分布位置)');xlabel('x');ylabel('y');zlabel('f(x,y)'); 
figure(3);
 max_global= max_global';
 index(:,1)=1:ECHO;
 plot(index(:,1),max_global(:,1),'b-')
 hold on;
 title('最优函数值变化趋势');
 xlabel('iteration');
 ylabel('f(x)');
 grid on;
 [c_max,i_max]=max(T0);
 maxpiont=[X(i_max,1),X(i_max,2)];
 maxvalue=cos(2*pi.*X(i_max,1)).*cos(2*pi.*X(i_max,2)).*exp(-((X(i_max,1).^2+X(i_max,2).^2)/10));
 

     
             
             
             
        
        
        
        