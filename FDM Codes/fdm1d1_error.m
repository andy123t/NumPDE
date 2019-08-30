% fdm1d1_error.m
% finite difference method for 1D problem
% -u''+u'=pi^2*sin(pi*x)+pi*cos(pi*x)  in [0,1]
% u(0)=0, u(1)=0 ;
% exact solution : u=sin(pi*x)
clear all
Nvec=[10 50 100 500 1000];    % Number of splits
Error=[];
for k=1:length(Nvec)
    N=Nvec(k);
    h=1/N;
    x=0:h:1;
    N=length(x)-1;
    A=diag((2/h^2)*ones(N-1,1))...
        +diag((1/(2*h)-1/h^2)*ones(N-2,1),1)...
        +diag((-1/(2*h)-1/h^2)*ones(N-2,1),-1);
    b=pi^2*sin(pi*x(2:N))+pi*cos(pi*x(2:N));
    u=A\b';
    u=[0;u;0];
    ue=sin(pi*x');
    error=max(abs(u-ue));
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1.5)
hold on,
plot(log10(Nvec), log10(Nvec.^(-2)), '--')
grid on,
xlabel('log_{10}N','fontsize', 16), ylabel('log_{10}Error','fontsize',16),
title('Convergence of Finite Difference Method','fontsize',14)
set(gca,'fontsize',14)
for i=1:length(Nvec)-1     % computating convergence order
    order(i)=-log(Error(i)/Error(i+1))/(log(Nvec(i)/Nvec(i+1)));
end
Error   
order
