clear all
close all


el = 0.0013 * 10 ^ -12; 
es = 10 * 10 ^ -12; 
Ee = 50 * 10 ^ -9; 
Ebf = 5 * 10 ^ -9;
n = 1000; d = 110; l = 4000; M = 100;
[m,d]= meshgrid(1:1:5,150:10:250);
f1= -Ebf .* m + el .* (d .^ 4) - (2 * m - 1) * Ee;
k=((n / (2 * 3.14))^(0.5))*((es ./ f1) .^ (0.5)) * M;
mesh(m, d, k);
xlabel('Head-Set Size', 'FontSize',14)
ylabel('Distance', 'FontSize',14)
zlabel('# of Clusters','FontSize',14)