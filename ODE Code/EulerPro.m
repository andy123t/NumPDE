% EulerPro.m
% Modified Euler method for the ODE model
% u'(x)=x^2+x-u, x in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(x)=-exp(-x)+x^2-x+1.
clear all;  clf
h=0.1;
x=0:h:1;                     % interval division
N=length(x)-1;
u(1)=0;                      % initial value
fun=@(x,u) x.^2+x-u;         % RHS
for n=1:N
    k1=fun(x(n),u(n));
    k2=fun(x(n+1),u(n)+h*k1);
    u(n+1)=u(n)+(h/2)*(k1+k2);
end
ue=-exp(-x)+x.^2-x+1;        % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
%title('Modified Euler Method','fontsize',12)
set(gca,'fontsize',12)
xlabel('x','fontsize', 16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  EulerPro.png
