% Trapezoidal.m
% Trapezoidal rule for the ODE model
% u'(t)=t^2+t-u, t in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(t)=-exp(-t)+t^2-t+1.
clear all; close all;
h=0.1;
t=0:h:1;                       % interval partition
N=length(t)-1;
un(1)=0;                       % initial value
for n=1:N
    un(n+1)=(2-h)/(2+h)*un(n)+h/(2+h)*(t(n)^2+t(n)+t(n+1)^2+t(n+1));
end
ue=-exp(-t)+t.^2-t+1;          % exact solution
plot(t,ue,'b-',t,un,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
%title('Trapezoidal rule','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% print -dpng -r600 Trapezoidal.png
% print -depsc2 Trapezoidal.eps
