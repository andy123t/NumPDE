% Trapezoidal.m
% Trapezoidal rule for the ODE model
% u'(x)=x^2+x-u, x in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(x)=-exp(-x)+x^2-x+1.
clear all;  clf
h=0.1;
x=0:h:1;                       % interval division
N=length(x)-1;
u(1)=0;                        % initial value
for n=1:N
    u(n+1)=(2-h)/(2+h).*u(n)+h/(2+h).*(x(n)^2+x(n)+x(n+1)^2+x(n+1));
end
ue=-exp(-x)+x.^2-x+1;          % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
%title('Trapezoidal rule','fontsize',12)
set(gca,'fontsize',12)
xlabel('x','fontsize',16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  Trapezoidal.png
