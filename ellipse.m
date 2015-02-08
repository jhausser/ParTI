function [] = ellipse(xc,yc,a,b,Q,style) 
%(xc,yc) Center of ellip 
% (a,b) Major and minor radius 
% Q - rotation matrix Orientation (degree) 
num = 20;      % #points -> smoothness 

theta = linspace(0,2*pi,num);
p(1,:) = a*cos(theta);
p(2,:) = b*sin(theta);

p = Q*p;

p(1,:) = p(1,:) + xc;
p(2,:) = p(2,:) + yc;

plot(p(1,:),p(2,:),style,'LineWidth',2)

end
