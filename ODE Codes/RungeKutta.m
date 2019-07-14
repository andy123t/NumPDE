% RungeKutta.m
% Runge-Kutta method for the ODE model
% u'=t^2+t-u,  t \in [0,1]
% Initial value : u(0)=0
% Exact : u(t)=-exp(-t)+t^2-t+1.
clear all
h=0.1;
x=0:h:1;                     % function interval
n=length(x)-1;
u(1)=0;                      % initial value
fun=@(t,u) t.^2+t-u;         % RHS
for i=1:n
    k1=fun(x(i),u(i));
    k2=fun(x(i)+h./2,u(i)+h.*k1/2);
    k3=fun(x(i)+h./2,u(i)+h.*k2/2);
    k4=fun(x(i)+h,u(i)+h.*k3);
    u(i+1)=u(i)+h.*(k1+2.*k2+2.*k3+k4)./6;
end
ue=-exp(-x)+x.^2-x+1;        % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1.5)
xlabel('x','fontsize', 16), ylabel('y','fontsize',16,'Rotation',0)
legend('Exact','Numerical','location','North')
title('Runge-Kutta Method','fontsize',14)
set(gca,'fontsize',14)
