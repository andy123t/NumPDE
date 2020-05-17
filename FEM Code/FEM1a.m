% FEM1a.m
% Finite Element Method
% -u_xx=-exp(x) in [0,1] 
% u(0)=0, u'(1)=e;
% exact solution: u=exp(x)-1.
clear all;  clf
N=10;
h=1/N;
X=0:h:1;
b=zeros(N,1);
A=zeros(N,N);
for i=2:N
    F1=@(x) N*(x-X(i-1)).*(-exp(x));     % handle
    R1=integral(F1,X(i-1),X(i));         % numerical integration
    F2=@(x) N*(X(i+1)-x).*(-exp(x));
    R2=integral(F2,X(i),X(i+1));
    b(i-1)=R1+R2;
end
F3=@(x) N*(x-X(N)).*(-exp(x));
b(N)=integral(F3,X(N),X(N+1))+exp(1); 
% Assemble matrix
a11=N;
a12=-N;
for i=1:N-1
    A(i,i)=2*a11;
    A(i,i+1)=a12;
    A(i+1,i)=a12;
end
A(N,N)=a11;
% solving the linear system
c=A\b;
un=vertcat(0,c); 
ue=exp(X)-1;
ue=ue';
Error=norm(un-ue)
% plot the figure
plot(X,ue,'b-',X,un,'r+','LineWidth',1)
%title('Numerical solutions vs Exact solutions','fontsize',14),
legend('Exact solutions','Numerical solutions','location','NorthWest')
set(gca,'fontsize',12)
xlabel('x','fontsize', 14), ylabel('u','fontsize',14,'Rotation',0)

% print -dpng -r600  FEM1a.png
