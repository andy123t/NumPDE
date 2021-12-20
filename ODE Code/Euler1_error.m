% Euler1_error.m
% Euler method for the ODE model
% u'(x)=x^2+x-u, x in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(x)=-exp(-x)+x^2-x+1.
clear all; close all;
Nvec=[10 50 100 500 1000];        % Number of partitions
Error=[];
fun=@(x,u) x.^2+x-u;              % RHS
for k=1:length(Nvec)
    N=Nvec(k);
    h=1/N;
    x=0:h:1;                      % interval division
    u(1)=0;                       % initial value
    for n=1:N
        u(n+1)=u(n)+h.*fun(x(n),u(n));
    end
    ue=-exp(-x)+x.^2-x+1;         % exact solution
    error=max(abs(u-ue));
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
%loglog(Nvec,Error,'ro-','LineWidth',1.5)
hold on
%loglog(Nvec, Nvec.^(-1), '--')
plot(log10(Nvec), log10(Nvec.^(-1)), '--')
grid on
%title('Convergence of Euler method','fontsize',12)
set(gca,'fontsize',12)
xlabel('log_{10}N','fontsize',14), ylabel('log_{10}Error','fontsize',14)

% add annotation of slope
ax = [0.57 0.53];
ay = [0.68 0.63];
annotation('textarrow',ax,ay,'String','slope = -1 ','fontsize',14)

% computating convergence order
for n=1:length(Nvec)-1
    order(n)=-log(Error(n)/Error(n+1))/(log(Nvec(n)/Nvec(n+1)));
end
Error
order

% print -dpng -r600 Euler1_error.png
% print -depsc2 Euler1_error.eps
