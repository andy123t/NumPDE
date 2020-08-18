% IRK2s2_error.m
% Implicit Runge-Kutta(Gauss method) 2 stage and order 4
% u'=u^2/2-1/2 in [0,1] with initial condition u(0)=0
% exact solution: ue=(1-exp(x))/(1+exp(x))
clear all;  close all;  clf
Nvec=[10 50 100 200 500];
Error=[];
for n=1:length(Nvec)
    N=Nvec(n);
    h=1/N;
    x=[0:h:1];
    u(1)=0;
    Y=[1;1];
    % Newton iteration
    for i=1:N
        k=u(i);
        tol=1;
        while tol>1.0e-6
            X=Y;
            D=[1-h*0.25*X(1), -h*(0.25-(sqrt(3))/6)*X(2);...
                -h*(0.25+(sqrt(3))/6)*X(1), 1-h*0.25*X(2)];
            F=[X(1)-k-h*(0.25*(0.5*(X(1))^2-0.5)+(0.25-(sqrt(3))/6)*(0.5*(X(2))^2-0.5));...
                X(2)-k-h*((0.25+(sqrt(3))/6)*(0.5*(X(1))^2-0.5)+0.25*(0.5*(X(2))^2-0.5))];
            Y=X-D\F;
            tol=norm(Y-X);
        end
        u(i+1)=k+(h/2)*(0.5*Y(1)^2-0.5+0.5*Y(2)^2-0.5);
    end
    ue=(1-exp(x))./(1+exp(x));  % exact solution
    error=max(abs(u-ue));       % maximum error
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
hold on,
plot(log10(Nvec), log10(Nvec.^(-4)), '--')
grid on,
%title('Convergence order of Gauss method ','fontsize',12)
set(gca,'fontsize',12)
xlabel('log_{10}N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)

% add annotation of slope
ax = [0.62 0.58];
ay = [0.72 0.66];
annotation('textarrow',ax,ay,'String','slope = -4 ','fontsize',14)

% computating convergence order
for i=1:length(Nvec)-1     % computating convergence order
    order(i)=-log(Error(i)/Error(i+1))/(log(Nvec(i)/Nvec(i+1)));
end
Error
order

% print -dpng -r600  IRK2s2_error.png
