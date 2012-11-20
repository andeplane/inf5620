function plotFE
close all;

phi1 = @(X) 0.5*(1-X);
phi2 = @(X) 0.5*(1+X);
X = linspace(-1,1,1000);
%h = 0.5;
h = pi/2;

%c = [h^2/6
%    12*(-35*h^3/72 + 7*h^2/12)/(7*h)
%    7*(-23*h^3/21 + 4*h^2/7)/(2*h)];
c = [6*(-4/pi + 4/3)/pi
    24*(-7/6 + 7/pi)/(7*pi)
    7*(-24/(7*pi) + 8/7)/pi];

x = linspace(0,pi,1000);
%f = x.*(1-x);
f = sin(x);

x1 = linspace(0,h,1000);
x2 = linspace(h,2*h,1000);

u1 = c(1) * phi1(X) + c(2) * phi2(X);
u2 = c(2) * phi1(X) + c(3) * phi2(X);

plot(x1,u1);
hold on
plot(x2,u2);
plot(x,f,'r');

b1 = phi1(X);
b2 = phi2(X);
c = b1.*b2;
figure
plot(c)


end