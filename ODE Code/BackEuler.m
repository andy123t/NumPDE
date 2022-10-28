% BackEuler.m
% Backward Euler Method for Nonlinear ODE:
% u'(t)=u-2t/u, t in [0,1]
% Initial condition: u(0)=1;
% Exact solution: u(t)=sqrt(2*t+1)
clear all; close all;
h=0.01;
t=0:h:1;
N=length(t)-1;
un(1)=1;
NI(1,N)=0;  % Record the number of iterations
for n=1:N
    % Newton iteration
    Xn=un(n);
    Xp=Xn;
    Xprev=0;
    while abs(Xp-Xprev) > abs(Xp)*1e-8
        Xprev=Xp;
        Xp=Xp-(Xp-h*Xp+2*h*t(n+1)/Xp-un(n))./(1-h-2*h*t(n+1)/(Xp^2));
        NI(n)=NI(n)+1;
    end
    un(n+1)=Xp;
end
ue=sqrt(2*t+1);    % exact solution
plot(t,ue,'b-',t,un,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
% title('Backward Euler method','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% computing error
error=max(abs(un-ue))

% print -dpng -r600 BackEuler.png
% print -depsc2 BackEuler.eps
