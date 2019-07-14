% Euler1_order.m
% Euler method for a first-order ODE
% u'(t)=t^2+t-u, t \in [0,2] 
% Initial value : u(0)=0 ;
% Exact solution : u(t)=-exp(-t)+t^2-t+1.
clear all
Nvec=[10 50 100 500 1000];  % Number of splits
Error=[];
fun=@(t,u) t.^2+t-u;            % RHS
for k=1:length(Nvec)
    N=Nvec(k);
    h=1/N;
    x=0:h:1;                                % function interval
    u(1)=0;                                 % initial value
    for i=1:N
        u(i+1)=u(i)+h.*fun(x(i),u(i));
    end
    ue=-exp(-x)+x.^2-x+1;        % exact solution
    error=max(abs(u-ue));
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','LineWidth',1.5)
%loglog(Nvec,Error,'ro-','LineWidth',1.5)
hold on,
%loglog(Nvec, Nvec.^(-1), '--')
plot(log10(Nvec), log10(Nvec.^(-1)), '--')
grid on,
xlabel('N','fontsize', 16), ylabel('Error','fontsize',16)
title('Convergence of Euler method','fontsize',14)
set(gca,'fontsize',14)
for i=1:length(Nvec)-1     % computating convergence order
    order(i)=-log(Error(i)/Error(i+1))/(log(Nvec(i)/Nvec(i+1)));
end
Error
order