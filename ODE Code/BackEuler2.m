% BackEuler2.m
% Backward Euler Method for Nonlinear ODE:
% u'(t)=-100*u^3+sin(u) in [0,1]
% Initial condition: u(0)=1.
clear all
h=0.001;
x=0:h:1;
N=length(x)-1;
u(1)=1;
NI(N)=0;  % Record the number of iterations
for n=1:N
    % Newton iteration
    Xn=u(n);
    Xp=Xn;
    Xprev=0;
    while abs(Xp-Xprev) > 1.0e-8*abs(Xp)
        Xprev=Xp;
        Xp=Xp-(Xp+100*h*Xp^3-h*sin(Xp)-u(n))./(1+300*h*Xp^2-h*cos(Xp));
        NI(n)=NI(n)+1;
    end
    u(n+1)=Xp;
end
plot(x,u,'r-.', 'LineWidth',1)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  BackEuler2.png
