% BackEuler2.m
% Backward Euler Method for Nonlinear ODE:
% u'(t)=-100*u^3+sin(u), t in [0,1]
% Initial condition: u(0)=1.
clear all; close all;
h=0.001;
t=0:h:1;
N=length(t)-1;
un(1)=1;
NI(N)=0;  % Record the number of iterations
for n=1:N
    % Newton iteration
    Xn=un(n);
    Xp=Xn;
    Xprev=0;
    while abs(Xp-Xprev) > abs(Xp)*1e-8
        Xprev=Xp;
        Xp=Xp-(Xp+100*h*Xp^3-h*sin(Xp)-un(n))./(1+300*h*Xp^2-h*cos(Xp));
        NI(n)=NI(n)+1;
    end
    un(n+1)=Xp;
end
plot(t,un,'r-.','LineWidth',1)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% print -dpng -r600 BackEuler2.png
% print -depsc2 BackEuler2.eps
