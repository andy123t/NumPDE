% ModiEuler.m
% Modified Euler method for the ODE model
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
    k2=fun(t(n+1),un(n)+h*k1);
    un(n+1)=un(n)+h/2*(k1+k2);
end
ue=-exp(-t)+t.^2-t+1;        % exact solution
plot(t,ue,'b-',t,un,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
%title('Modified Euler Method','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% print -dpng -r600 ModiEuler.png
% print -depsc2 ModiEuler.eps
