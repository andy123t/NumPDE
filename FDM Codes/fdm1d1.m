% fdm1d1.m
% finite difference method for 1D problem
% -u''+u'=pi^2*sin(pi*x)+pi*cos(pi*x)  in [0,1]
% u(0)=0, u(1)=0 ;
% exact solution : u=sin(pi*x)
clear all
h=0.05;
x=0:h:1;
N=length(x)-1;
A=diag((2/h^2)*ones(N-1,1))...
    +diag((1/(2*h)-1/h^2)*ones(N-2,1),1)...
    +diag((-1/(2*h)-1/h^2)*ones(N-2,1),-1);
b=pi^2*sin(pi*x(2:N))+pi*cos(pi*x(2:N));
u=A\b';
u=[0;u;0];
ue=sin(pi*x)';
plot(x,ue,'b-',x,u,'r+','LineWidth',1.5)
Error=max(abs(u-ue))
xlabel('x','fontsize', 16), ylabel('y','fontsize',16,'Rotation',0)
legend('Exact ','Numerical','location','North')
title('Finite Difference Method','fontsize',14)
set(gca,'fontsize',14)
