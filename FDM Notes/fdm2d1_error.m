% fdm2d1_error.m
% finite difference method for 2D problem
% -d^2u/dx^2-d^2u/dy^2=f(x,y)
% f(x,y)=-2*pi^2*exp(pi*(x+y))*(sin(pi*x)*cos(pi*y)+cos(pi*x)*sin(pi*y))
% exact solution: ue=exp(pi*x+pi*y)*sin(pi*x)*sin(pi*y)
clear all
Nvec=2.^[3:10];
Err=[];
for n=Nvec
   h=1/n;
   x=[0:h:1]';
   y=[0:h:1]';
   N=length(x)-1;
   M=length(y)-1;
   [X,Y]=meshgrid(x,y);
   X=X(2:M,2:N);
   Y=Y(2:M,2:N);
   % generate the matrix of RHS
   f=-2*pi^2*exp(pi*X+pi*Y).*(sin(pi*X).*cos(pi*Y)+cos(pi*X).*sin(pi*Y));
   % constructing the coefficient matrix
   e=ones(N-1,1);
   C=1/h^2*spdiags([-e 4*e -e],[-1 0 1],N-1,N-1);
   D=-1/h^2*eye(N-1);
   e=ones(M-1,1);
   A=kron(eye(M-1),C)+kron(spdiags([e e],[-1 1],M-1,M-1),D);
   % solving the linear system
   f=f';
   u=zeros(M+1,N+1);
   u(2:M,2:N)=reshape(A\f(:),N-1,M-1)';
   u(:,1)=0;
   u(:,end)=0;
   ue=zeros(M+1,N+1);
   % numerical solution
   ue(2:M,2:N)=exp(pi*X+pi*Y).*sin(pi*X).*sin(pi*Y);
   err=max(max(abs(u-ue)));     % maximum error
   Err=[Err,err];
end
plot(log10(Nvec),log10(Err),'ro-','MarkerFaceColor','w','LineWidth',1.5)
hold on,
plot(log10(Nvec), log10(Nvec.^(-2)), '--')
grid on,
xlabel('log_{10}N','fontsize', 16), ylabel('log_{10}Error','fontsize',16),
title('Convergence of Finite Difference Method','fontsize',14)
set(gca,'fontsize',14)

for i=1:length(Nvec)-1     % computating convergence order
   order(i)=-log(Err(i)/Err(i+1))/(log(Nvec(i)/Nvec(i+1)));
end
Err
order
