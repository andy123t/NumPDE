% fdm1d2.m
% finite difference method for 1D problem
% -(xu')'+x*u'=pi^2*x*sin(pi*x)-pi*cos(pi*x)+pi*x*cos(pi*x) in [0,1]
% u(0)=0, u(1)=0 ;
% exact solution : u=sin(pi*x)
clear all;  clf
h=0.05;
x=0:h:1;
N=length(x)-1;
A=diag(2*x(2:N)./h^2)+diag(x(2:N-1)./(2*h)-(x(2:N-1)+0.5*h)./h^2,1)...
    +diag(-x(3:N)./(2*h)-(x(3:N)-0.5*h)./h^2,-1);
b=pi^2*x(2:N).*sin(pi*x(2:N))+pi*(x(2:N)-1).*cos(pi*x(2:N));
u=A\b';
u=[0;u;0];
ue=sin(pi*x');
plot(x,ue,'b-',x,u,'r+','LineWidth',1)
Error=max(abs(u-ue))
legend('Exact ','Numerical','location','NorthEast')
%title('Finite Difference Method','fontsize',12)
set(gca,'fontsize',12)
xlabel('x','fontsize', 16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  fdm1d2.png
