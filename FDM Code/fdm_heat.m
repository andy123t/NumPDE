% fdm_heat.m
% forward difference scheme for heat equation
% u_t=u_{xx}, (x,t) in (0,1)x(0,1],
% u(x,0)=exp(x), x in [0,1],
% u(0,t)=exp(t), u(1,t)=exp(1+t), t in (0,1]
% exact solution: u(x,t)=exp(x+t)
clear all; close all;
a=1;
h=0.05; x=0:h:1;
tau=0.00125; t=0:tau:1;
r=a*tau/h^2;
M=length(x)-1; N=length(t)-1;
% constructing the coefficient matrix
e=r*ones(M-1,1);
A=spdiags([e 1-2*e e],[-1 0 1],M-1,M-1);
% setting initial and boundary conditions
un=zeros(M+1,N+1);
un(:,1)=exp(x);
un(1,:)=exp(t); 
un(end,:)=exp(1+t);
for n=1:N
    un(2:M,n+1)=A*un(2:M,n);
    un(2,n+1)=un(2,n+1)+r*un(1,n);
    un(M,n+1)=un(M,n+1)+r*un(end,n);
end
% plot the figure
mesh(t(1:20:end),x,un(:,1:20:end))
set(gca,'fontsize',12)
xlabel('t','fontsize', 14) 
ylabel('x','fontsize',14)
zlabel('u','fontsize',14)

% calculating maximum error
[T,X]=meshgrid(t,x);
ue=exp(X+T);
Error=max(max(abs(ue-un)))

% print -dpng -r600 fdm_heat.png
% print -depsc2 fdm_heat.eps
