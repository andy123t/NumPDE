% fdm_heat.m
% forward difference scheme for heat equation
% u_t=u_{xx}, (x,t) in (0,1)x(0,1],
% u(x,0)=exp(x), x in [0,1],
% u(0,t)=exp(t), u(1,t)=exp(1+t), t in (0,1]
% exact solution: u(x,t)=exp(x+t)
clear all; close all;
a=1;
h=0.05; x=[0:h:1];
tau=0.00125; t=[0:tau:1];
r=a*tau/h^2;
M=length(x)-1; N=length(t)-1;
% constructing the coefficient matrix
e=r*ones(M-1,1);
A=spdiags([e 1-2*e e],[-1 0 1],M-1,M-1);
% setting initial and boundary conditions
u=zeros(M+1,N+1);
u(:,1)=exp(x);
u(1,:)=exp(t); 
u(end,:)=exp(1+t);
for n=1:N
    u(2:M,n+1)=A*u(2:M,n);
    u(2,n+1)=u(2,n+1)+r*u(1,n);
    u(M,n+1)=u(M,n+1)+r*u(end,n);
end
% plot the figure
mesh(t(1:20:end),x,u(:,1:20:end))
set(gca,'fontsize',12)
xlabel('t','fontsize', 14) 
ylabel('x','fontsize',14)
zlabel('u','fontsize',14,'Rotation',0)

% calculating maximum error
[T X]=meshgrid(t,x);
ue=exp(X+T);
Error=max(max(abs(ue-u)))

% print -dpng -r600  fdm_heat.png
