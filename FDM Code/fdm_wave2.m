% fdm_wave2.m
% finite difference method for wave equation
% u_{tt}=u_{xx}, (x,t) in (0,1)x(0,T],
% u(x,0)=0, u_t(x,0)=x, x in [0,1],
% u(0,t)=0, u(1,t)=sin(t), t in (0,1].
% exact solution: u(x,t)=sin(4*pi*x)*cos(4*pi*t)
%		+sin(8*pi*x)*sin(8*pi*t)/(8*pi);
clear all; close all;
a=1; 
T=1;   % time
h=0.0025; x=0:h:1;
tau=0.0025; t=0:tau:T;
r=a*tau/h; 
M=length(x)-1; N=length(t)-1;
[T,X]=meshgrid(t,x);
% constructing the coefficient matrix
e=r^2*ones(M-1,1);
A=spdiags([e 2*(1-e) e],[-1 0 1],M-1,M-1);
% setting initial and boundary conditions
un=zeros(M+1,N+1);
un(:,1)=sin(4*pi*x');
un(:,2)=sin(4*pi*x')+tau*sin(8*pi*x');
un(1,:)=0; un(end,:)=0;
for n=2:N
    un(2:M,n+1)=A*un(2:M,n)-un(2:M,n-1);
    un(2,n+1)=un(2,n+1)+r^2*un(1,n);
    un(M,n+1)=un(M,n+1)+r^2*un(end,n);
end
% plot the figure
mesh(t,x,un), view(20,40)
set(gca,'fontsize',12)
xlabel('t','fontsize', 14) 
ylabel('x','fontsize',14)
zlabel('u','fontsize',14)

% calculating maximum error
ue=sin(4*pi*X).*cos(4*pi*T)+sin(8*pi*X).*sin(8*pi*T)/(8*pi);

disp('t = 1, 2, 3, 4, 5')
Error=max(abs(ue(:,1/tau+1:1/tau:end)-un(:,1/tau+1:1/tau:end)))

% print -dpng -r600 fdm_wave2.png
% print -depsc2 fdm_wave2.eps
