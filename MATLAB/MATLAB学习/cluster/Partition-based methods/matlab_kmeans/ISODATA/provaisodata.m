clear all
close all
clc
load('\dades.mat')
ON=15;  %Umbral del n�mero de elementos para la eliminaci�n de un agrupamiento. 
OC=10;  %Umbral de distancia para la uni�n de agrupamientos.
OS=7;  %Umbral de desviaci�n t�pica para la divisi�n de un agrupamiento.
k=4;   %N�mero (m�ximo) de agrupamientos.
L=2;   %M�ximo n�mero de agrupamientos que pueden mezclarse en una sola iteraci�n.
I=10;  %M�ximo n�mero de iteraciones permitidas.
NO=1;  %Parametro extra para responder automaticamente que no a la peticion de cambial algun parametro.
min=50; %Minima distancia que un punto debe de estar de cada centro. Si no deseas eliminar ningun punto
        % dale un valor elevado.
%%%%%%%%%%%%%%%%%%%%%
%  Funcion ISODATA  %
%%%%%%%%%%%%%%%%%%%%%
[centro, Xcluster, Ycluster, A, clustering]=isodata(X, Y, k, L, I, ON, OC, OS, NO, min);
clc;
fprintf('Numero de agrupaciones: %d',A);

% Presentacion de resultados por pantalla.

% Creamos los colores.
colr=zeros(A,3);
for i=1:A
    colr(i,:)=rand(1,3);
end;

% Representamos la informacion.
figure;
hold on;
for i=1:A,
    n=find(clustering==i);
    p=plot(X(n), Y(n),'.');set(p,'Color',colr(i,:));title(A)
end;

%plot(centro(:,1), centro(:,2), 'g.');

clc;
fprintf('Numero de agrupaciones: %d',A);
% Borramos variables temporales.
clear n;clear i;clear p;clear colr;clear NO;