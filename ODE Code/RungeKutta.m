% RungeKutta.m
% Runge-Kutta method for the ODE model
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
    k2=fun(x(n)+h./2,u(n)+h.*k1/2);
    k3=fun(x(n)+h./2,u(n)+h.*k2/2);
    k4=fun(x(n)+h,u(n)+h.*k3);
    u(n+1)=u(n)+h.*(k1+2.*k2+2.*k3+k4)./6;
end
ue=-exp(-x)+x.^2-x+1;        % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
% title('Runge-Kutta Method','fontsize',12)
set(gca,'fontsize',12)
xlabel('x','fontsize',16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  RungeKutta.png
