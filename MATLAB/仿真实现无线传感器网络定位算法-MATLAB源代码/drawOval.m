%| function [h, j] = drawOval(centerX, centerY, R, lineStyle, markChar, markSize, displayOn, sigma)
%|
%| INPUTS:
%|   (centerX, centerY) define the center of the oval.
%|   R is a covariance matrix of which we want to plot the 1-sigma
%|     uncertainty ellipse.
%|   The oval is drawn with axis lengths equal to the std deviation on each axis,
%|     and axis equal to the eigenvectors of R. 
%|   The 'linestyle' string must have the first character as a color (eg 'r')
%|
%|  Neal Patwari, 10/1/03  (http://www.engin.umich.edu/~npatwari/)
%|

function [h, j] = drawOval(centerX, centerY, R, lineStyle, markChar, markSize, displayOn, sigma)

if ~exist('sigma')
    sigma = 1;
end
if ~exist('displayOn')
   displayOn = 0;
end

if abs(R(1,2)) < eps,
   alpha1 = 0;
   alpha2 = pi/2;
   rads   = sqrt(diag(R));
   v      = [1 0; 0 1];
elseif max(max(R))==inf || max(max(R))== NaN,
    disp('Unable to display oval for R')
    disp(R);
    alpha1=0;
    alpha2=pi/2;
    rads=[eps eps];
else
   [v, d] = eig(R);
   alpha1 = atan2(v(2,1),v(1,1));
   alpha2 = atan2(v(2,2),v(1,2));
   rads   = sqrt(diag(d));
end
deltaT = 2*pi/100;
theta  = 0:deltaT:2*pi;
x      = sigma*rads(1).*cos(theta);
y      = sigma*rads(2).*sin(theta);
xprime = centerX + x.*cos(alpha1) - y.*sin(alpha1);
yprime = centerY + x.*sin(alpha1) + y.*cos(alpha1);
h = plot(xprime, yprime, lineStyle);
set(h,'LineWidth',2)
j = plot(centerX, centerY, sprintf('%s%s',lineStyle(1), markChar));
if exist('markSize'),
   set(j,'MarkerSize',markSize)
   set(j,'MarkerFaceColor',lineStyle(1))
end
if displayOn,
    R
    v
    rads
    disp([alpha1, alpha2])
   axis1 = [v(:,1), -v(:,1)] * rads(1) * 1.1;
   axis2 = [v(:,2), -v(:,2)] * rads(2) * 1.1;
   plot(centerX + axis1(1,:), centerY + axis1(2,:), 'k:')
   plot(centerX + axis2(1,:), centerY + axis2(2,:), 'k:')
   set(gca,'ylim',[-0.1 1.1])
   set(gca,'xlim',[-0.1 1.1])
   grid
   set(gca,'DataAspectRatio',[1,1,1])
end

