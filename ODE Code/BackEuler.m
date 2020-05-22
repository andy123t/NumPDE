% BackEuler.m
% Backward Euler Method for Nonlinear ODE:
% u'(x)=u-2x/u in [0,1]
% Initial condition: u(0)=1;
% Exact solution: u=sqrt(2*x+1)
clear all
h=0.05;
x=0:h:1;
N=length(x)-1;
u(1)=1;
NI(N)=0;           % Record the number of iterations
for n=1:N
    % Newton iteration
    X(1)=u(n); tol=1.0e-6;
    for k=1:100
        X(k+1)=X(k)-(X(k)-h*X(k)+2*h*x(n+1)/X(k)-u(n))./(1-h-2*h*x(n+1)/(X(k)^2));
        NI(n)=NI(n)+1;
        if abs(X(k+1)-X(k))<tol
            u(n+1)=X(k+1);
            break
        end
    end
end
ue=sqrt(2*x+1);    % exact solution
plot(x,ue,'b-',x,u,'r+','LineWidth',1)
legend('Exact ','Numerical','location','North')
% title('Backward Euler method','fontsize',12)
set(gca,'fontsize',12)
xlabel('t','fontsize', 16), ylabel('u','fontsize',16,'Rotation',0)
% computing error
error=max(abs(u-ue))