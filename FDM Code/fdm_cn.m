% fdm_cn.m
% Crank-Nicolson scheme for heat equation
% u_t=u_{xx}+f(x,t), (x,t) in (0,1)x(0,1],
% u(x,0)=sin(x), x in [0,1],
% u(0,t)=sin(t), u(1,t)=sin(1+t), t in (0,1]
% f(x,t)=cos(x+t)+sin(x+t),
% exact solution: u(x,t)=sin(x+t)
clear all; close all;
a=1;
h=0.05; x=[0:h:1];
tau=0.001; t=[0:tau:1];
r=a*tau/h^2;
M=length(x)-1; N=length(t)-1;
% constructing the coefficient matrix
e=r*ones(M-1,1);
A=spdiags([-e 2+2*e -e],[-1 0 1],M-1,M-1);
B=spdiags([e 2-2*e e],[-1 0 1],M-1,M-1);
% setting initial and boundary conditions
u=zeros(M+1,N+1);
u(:,1)=sin(x);
u(1,:)=sin(t); 
u(end,:)=sin(1+t);
for n=1:N
    F=tau*cos(x(2:M)'+t(n))+tau*sin(x(2:M)'+t(n))...
        +tau*cos(x(2:M)'+t(n+1))+tau*sin(x(2:M)'+t(n+1));
    F(1)=F(1)+r*u(1,n)+r*u(1,n+1);
    F(M-1)=F(M-1)+r*u(end,n)+r*u(end,n+1);
    % solving the system
    u(2:M,n+1)=A\B*u(2:M,n)+A\F;
end
% plot the figure
mesh(t(1:20:end),x,u(:,1:20:end))
set(gca,'fontsize',12)
xlabel('t','fontsize', 14) 
ylabel('x','fontsize',14)
zlabel('u','fontsize',14,'Rotation',0)

% calculating maximum error
[T X]=meshgrid(t,x);
ue=sin(X+T);

Error=max(max(abs(ue-u)))

% print -dpng -r600  fdm_cn.png
