% BDF2_error.m
% BDF method for the ODE model
% u'(t)=t^2+t-u, t in [0,1] 
% Initial condition: u(0)=0 ;
% Exact solution: u(t)=-exp(-t)+t^2-t+1.
clear all; close all;
Nvec=[10 50 100 500 1000];   % Number of partitions
Error=[];
for k=1:length(Nvec)
    N=Nvec(k);
    h=1/N;
    t=0:h:1;                 % interval partition
    un(1)=0;                 % initial value
    % Trapezoidal rule to compute un(2)
    un(2)=2/(2+h)*un(1)+h/(2+h)*(t(1)^2+t(1)+un(1))+h/(2+h)*(t(2)^2+t(2));
    % BDF2 method
    for n=1:N-1
        un(n+2)=4/(3+2*h)*un(n+1)-1/(3+2*h)*un(n)+2*h/(3+2*h)*(t(n+2)^2+t(n+2));
    end
    ue=-exp(-t)+t.^2-t+1;    % exact solution
    error=max(abs(un-ue));
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
%loglog(Nvec,Error,'ro-','LineWidth',1.5)
hold on
%loglog(Nvec,Nvec.^(-2),'--')
plot(log10(Nvec),log10(Nvec.^(-2)),'--')
grid on
%title('Convergence of BDF2 method','fontsize',12)
set(gca,'fontsize',12)
xlabel('log_{10}N','fontsize',14), ylabel('log_{10}Error','fontsize',14)

% add annotation of slope
ax = [0.57 0.53];
ay = [0.68 0.63];
annotation('textarrow',ax,ay,'String','slope = -2','fontsize',14)

% computing convergence order
for n=1:length(Nvec)-1
    order(n)=-log(Error(n)/Error(n+1))/(log(Nvec(n)/Nvec(n+1)));
end
Error
order

% print -dpng -r600 BDF2_error.png
% print -depsc2 BDF2_error.eps
