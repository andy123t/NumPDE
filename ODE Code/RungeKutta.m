% RungeKutta.m
% Runge-Kutta method for the ODE model
% u'(t)=t^2+t-u, t in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(t)=-exp(-t)+t^2-t+1.
clear all; close all;
h=0.1;
t=0:h:1;                     % interval partition
N=length(t)-1;
un(1)=0;                     % initial value
fun=@(x,u) x.^2+x-u;         % RHS
for n=1:N
    k1=fun(t(n),un(n));
    k2=fun(t(n)+h/2,un(n)+h/2*k1);
    k3=fun(t(n)+h/2,un(n)+h/2*k2);
    k4=fun(t(n)+h,un(n)+h*k3);
    un(n+1)=un(n)+h/6*(k1+2*k2+2*k3+k4);
end
ue=-exp(-t)+t.^2-t+1;        % exact solution
plot(t,ue,'b-',t,un,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
% title('Runge-Kutta Method','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% print -dpng -r600 RungeKutta.png
% print -depsc2 RungeKutta.eps
