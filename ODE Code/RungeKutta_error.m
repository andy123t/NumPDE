% RungeKutta_error.m
% Runge-Kutta method for the ODE model
% u'(t)=t^2+t-u, t \in [0,2] 
% Initial value : u(0)=0 ;
% Exact solution : u(t)=-exp(-t)+t^2-t+1.
clear all;  clf
Nvec=[10 50 100 500 1000];        % Number of splits
Error=[];
fun=@(t,u) t.^2+t-u;              % RHS
for k=1:length(Nvec)
    N=Nvec(k);
    h=1/N;
    x=0:h:1;                      % function interval
    u(1)=0;                       % initial value
    for n=1:N
        k1=fun(x(n),u(n));
        k2=fun(x(n)+h./2,u(n)+h.*k1/2);
        k3=fun(x(n)+h./2,u(n)+h.*k2/2);
        k4=fun(x(n)+h,u(n)+h.*k3);
        u(n+1)=u(n)+h.*(k1+2.*k2+2.*k3+k4)./6;
    end
    ue=-exp(-x)+x.^2-x+1;         % exact solution
    error=max(abs(u-ue));
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(log10(Nvec), log10(Nvec.^(-4)), '--')
grid on
%title('Convergence of Runge-Kutta Method','fontsize',12)
set(gca,'fontsize',12)
xlabel('log_{10}N','fontsize', 16), ylabel('log_{10}Error','fontsize',16)

% add annotation of slope
ax = [0.57 0.53];
ay = [0.68 0.63];
annotation('textarrow',ax,ay,'String','slope = -4 ','fontsize',14)

% computating convergence order
for n=1:length(Nvec)-1
    order(n)=-log(Error(n)/Error(n+1))/(log(Nvec(n)/Nvec(n+1)));
end
Error
order

% print -dpng -r600  RungeKutta_error.png
