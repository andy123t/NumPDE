% fdm_wave.m
% finite difference method for wave equation
% u_{tt}=u_{xx}+f(x,t), (x,t) in (0,1)x(0,1],
% u(x,0)=0, u_t(x,0)=x, x in [0,1],
% u(0,t)=0, u(1,t)=sin(t), t in (0,1].
% f(x,t)=cos(x+t)+sin(x+t),
% exact solution: u(x,t)=sin(x+t)
clear all; close all;
a=1; 
h=0.05; x=[0:h:1];
tau=0.05; t=[0:tau:1];
r=a*tau/h; 
M=length(x)-1; N=length(t)-1;
[T X]=meshgrid(t,x);
% constructing the coefficient matrix
e=r^2*ones(M-1,1);
A=spdiags([e 2*(1-e) e],[-1 0 1],M-1,M-1);
% setting initial and boundary conditions
u=zeros(M+1,N+1);
u(:,1)=0; u(:,2)=tau*x;
u(1,:)=0; u(end,:)=sin(t);
for n=2:N
    u(2:M,n+1)=A*u(2:M,n)-u(2:M,n-1)+ ...
        tau^2*(T(2:M,n).^2-X(2:M,n).^2).*sin(X(2:M,n).*T(2:M,n));
    u(2,n+1)=u(2,n+1)+r^2*u(1,n);
    u(M,n+1)=u(M,n+1)+r^2*u(end,n);
end
% plot the figure
mesh(t,x,u), view(20,40)
set(gca,'fontsize',12)
xlabel('t','fontsize', 14) 
ylabel('x','fontsize',14)
zlabel('u','fontsize',14)

% calculating maximum error
ue=sin(X.*T);
Error=max(max(abs(ue-u)))

% print -dpng -r600 fdm_wave.png
% print -depsc2 fdm_wave.eps
