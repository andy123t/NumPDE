% BackEuler_error.m
% Backward Euler Method for Nolinear ODE:
% u'(x)=u-2x/u in [0,1]
% Initial condition: u(0)=1;
% Exact solution: u=sqrt(2*x+1)
clear all
Nvec=[10 20 100 200 1000];
%hv=[0.1 0.05 0.01 0.005 0.001];
Error=[];
for k=1:length(Nvec)
N=Nvec(k);
h=1/N;
x=0:h:1;
N=length(x)-1;
u(1)=1;
for n=1:N
    % Newton iteration
    X(1)=u(n); tol=1.0e-10;
    for k=1:100
        X(k+1)=X(k)-(X(k)-h*X(k)+2*h*x(n+1)/X(k)-u(n))./(1-h-2*h*x(n+1)/(X(k)^2));
        if abs(X(k+1)-X(k))<tol
            u(n+1)=X(k+1);
            break
        end
    end
end
ue=sqrt(2*x+1);
error=max(abs(u-ue));
Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(log10(Nvec), log10(Nvec.^(-1)), '--')
grid on
%title('Rate of Convergence','fontsize',12)
set(gca,'fontsize',14)
xlabel('log_{10}N','fontsize', 14),ylabel('log_{10}Error','fontsize', 14)

% add annotation of slope
ax = [0.57 0.53];
ay = [0.68 0.63];
annotation('textarrow',ax,ay,'String','slope = -1 ','fontsize',14)

% computating convergence order
for n=1:length(Nvec)-1     % computating convergence order
    order(n)=-log(Error(n)/Error(n+1))/(log(Nvec(n)/Nvec(n+1)));
end
Error
order

% print -dpng -r600  BackEuler_error.png
