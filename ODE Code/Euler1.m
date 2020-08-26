% Euler1.m
% Euler method for the ODE model
% u'(x)=x^2+x-u, x in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(x)=-exp(-x)+x^2-x+1.
clear all;  clf
h=0.1;
x=0:h:1;                       % interval division
N=length(x)-1;
u(1)=0;                        % initial value
fun=@(t,u) t.^2+t-u;           % RHS
for n=1:N
    u(n+1)=u(n)+h.*fun(x(n),u(n));
end
ue=-exp(-x)+x.^2-x+1;          % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
%title('Euler method','fontsize',12)
set(gca,'fontsize',12)
xlabel('x','fontsize', 16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  Euler1.png
