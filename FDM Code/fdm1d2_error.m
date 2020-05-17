% fdm1d2_error.m
% finite difference method for 1D problem
% -(xu')'+x*u'=pi^2*x*sin(pi*x)-pi*cos(pi*x)+pi*x*cos(pi*x) in [0,1]
% u(0)=0, u(1)=0 ;
% exact solution : u=sin(pi*x)
clear all;  clf
Nvec=[10 20 50 100 200 500 1000];    % Number of splits
Error=[];
for k=1:length(Nvec)
    N=Nvec(k);
    h=1/N;
    x=0:h:1;
    N=length(x)-1;
    A=diag(2*x(2:N)./h^2)+diag(x(2:N-1)./(2*h)-(x(2:N-1)+0.5*h)./h^2,1)...
        +diag(-x(3:N)./(2*h)-(x(3:N)-0.5*h)./h^2,-1);
    b=pi^2*x(2:N).*sin(pi*x(2:N))+pi*(x(2:N)-1).*cos(pi*x(2:N));
    u=A\b';
    u=[0;u;0];
    ue=sin(pi*x');
    error=max(abs(u-ue));
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(log10(Nvec), log10(Nvec.^(-2)), '--')
grid on
%title('Convergence of Finite Difference Method','fontsize',14)
set(gca,'fontsize',14)
xlabel('log_{10}N','fontsize', 14), ylabel('log_{10}Error','fontsize',14),

% add annotation of slope
ax = [0.57 0.53];
ay = [0.68 0.63];
annotation('textarrow',ax,ay,'String','slope = -2 ','fontsize',14)

% computating convergence order
for i=1:length(Nvec)-1
    order(i)=-log(Error(i)/Error(i+1))/(log(Nvec(i)/Nvec(i+1)));
end
Error
order

% print -dpng -r600  fdm1d2_error.png
