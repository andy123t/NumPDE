% BDF2.m
% BDF method for the ODE model
% u'(t)=t^2+t-u, t in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(t)=-exp(-t)+t^2-t+1.
clear all; close all;
h=0.1;
t=0:h:1;                       % interval partition
N=length(t)-1;
un=zeros(1,N+1);
un(1)=0;                       % initial value
% Implicit Euler method to compute un(2)
% un(2)=1/(1+h)*un(1)+h/(1+h)*(t(2)^2+t(2));
% Trapezoidal rule to compute un(2)
un(2)=2/(2+h)*un(1)+h/(2+h)*(t(1)^2+t(1)+un(1))+h/(2+h)*(t(2)^2+t(2));
% BDF2 method
for n=1:N-1
    un(n+2)=4/(3+2*h)*un(n+1)-1/(3+2*h)*un(n)+2*h/(3+2*h)*(t(n+2)^2+t(n+2));
end
ue=-exp(-t)+t.^2-t+1;          % exact solution
plot(t,ue,'b-',t,un,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
% title('BDF2 method','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% computing error
error=max(abs(un-ue))

% print -dpng -r600 BDF2.png
% print -depsc2 BDF2.eps
