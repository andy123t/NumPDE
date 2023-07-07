% BDF2Non.m
% BDF method for the nonlinear ODE model
% u'(t)=u-2t/u, t in [0,1]
% Initial condition: u(0)=1;
% Exact solution: u(t)=sqrt(2*t+1)
clear all; close all;
h=0.1;
t=0:h:1;          % interval partition
N=length(t)-1;
un=zeros(1,N+1);
un(1)=1;          % initial value
NI(1,N)=0;        % Record the number of iterations
% Trapezoidal rule to compute un(2)
X=un(1);
Xprev=0;
while abs(X-Xprev) > abs(X)*1e-14
    Xprev=X;
    X=X-(X-h/2*X+h*t(2)/X-un(1)-h/2*(un(1)...
        -2*t(1)/un(1)))./(1-h/2-h*t(2)/(X^2));
    un(2)=X;
    NI(1)=NI(1)+1;
end
% BDF2 method
for n=1:N-1
    X=un(n+1);
    Xprev=0;
    while abs(X-Xprev) > abs(X)*1e-14
    Xprev=X;
    X=X-(X-2/3*h*X+4/3*h*t(n+2)/X-4/3*un(n+1)...
        +1/3*un(n))./(1-2/3*h-4/3*h*t(n+2)/(X^2));
    NI(n+1)=NI(n+1)+1;
    end
    un(n+2)=X;
end
ue=sqrt(2*t+1);    % exact solution
plot(t,ue,'b-',t,un,'r+','LineWidth',1)
legend('Exact','Numerical','location','northwest')
% title('BDF2 method','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize',16), ylabel('u','fontsize',16)

% computing error
error=max(abs(un-ue))

% print -dpng -r600 BDF2Non.png
% print -depsc2 BDF2Non.eps
