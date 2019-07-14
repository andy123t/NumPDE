% Euler1.m
% Euler method for a first-order ODE
% u'(t)=t^2+t-u, t \in [0,1] 
% Initial value : u(0)=0 ;
% Exact solution : u(t)=-exp(-t)+t^2-t+1.
clear all
h=0.1;
x=0:h:1;                               % function interval
n=length(x)-1;
u(1)=0;                                % initial value
fun=@(t,u) t.^2+t-u;           % RHS
for i=1:n
    u(i+1)=u(i)+h.*fun(x(i),u(i));
end
ue=-exp(-x)+x.^2-x+1;        % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1.5)
xlabel('x','fontsize', 16), ylabel('y','fontsize',16,'Rotation',0)
legend('Exact ','Numerical','location','North')
title('Euler Method','fontsize',14)
set(gca,'fontsize',14)
