% BackEuler2.m
% Backward Euler Method for Nonlinear ODE:
% u'(t)=-100*u^3+sin(u) in [0,1]
% Initial condition: u(0)=1.
clear all
h=0.01;
x=0:h:1;
N=length(x)-1;
u(1)=1;
NI(N)=0;       % Record the number of iterations
for n=1:N
    % Newton iteration
    X(1)=u(n); tol=1.0e-8;
    for k=1:100
        X(k+1)=X(k)-(X(k)+100*h*X(k)^3-h*sin(X(k))-u(n))./(1+300*h*X(k)^2-h*cos(X(k)));
        NI(n)=NI(n)+1;
        if abs(X(k+1)-X(k))<tol
            u(n+1)=X(k+1);
            break
        end
    end
end
plot(x,u,'r-.', 'LineWidth',1)
set(gca,'fontsize',12)
xlabel('t','fontsize', 16), ylabel('u','fontsize',16,'Rotation',0)

% print -dpng -r600  BackEuler2.png
