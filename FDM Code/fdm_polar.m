% fdm_polar.m
% finite difference method for polar coordinate problem
% -{\Delta}_{r,theta}u(r,theta)=f(r,theta), [0,1]x[0,2*pi]
% u(1,theta)=0,  theta in [0,2*pi]
% exact solution: u=(1-r^2)*r*(sin(theta)+cos(theta))/4.
clear all; close all;
N=50;
M=100;
dr=1/(N+1/2);
% dr=0.02;
dthe=2*pi/M;
r=(dr/2:dr:1)';
the=(0:dthe:2*pi)';
% N=length(r)-1;
% M=length(the);
% generate coordinates on the grid
[The,R]=meshgrid(the,r);
% generate the matrix of RHS
The1=The(1:N,2:end); R1=R(1:N,2:end);
f=2*R1.*(sin(The1)+cos(The1));
f=f';
% constructing the coefficient matrix
e=[2/dr^2+8/(dr^2*dthe^2);(2./(dthe^2*r(2:N).^2)+2/dr^2)];
e1=[-2/dr^2;-(r(2:N-1)+dr/2)./(dr^2*r(2:N-1))];
e2=-(r(2:N)-dr/2)./(dr^2*r(2:N));
C=diag(e)+diag(e1,1)+diag(e2,-1);
D=diag([-4/(dr^2*dthe^2);-1./(dthe^2*r(2:N).^2)]);
E=diag([-4/(dr^2*dthe^2);-1./(dthe^2*r(2:N).^2)]);
A=kron(eye(M),C)+kron(diag(ones(M-1,1),1)+diag(1,1-M),D)...
    +kron(diag(ones(M-1,1),-1)+diag(1,M-1),E);
% solving the linear system
f=f';
% u=zeros(M+1,N+1);
uh=reshape(A\f(:),N,M);
un=[uh;zeros(1,M)];
un=[un(:,end),un];
ue=(1-R.^2).*R.*(sin(The)+cos(The))/4;
% error on the node of mesh
Err=abs(un(1:N,2:end)-ue(1:N,2:end));
% compute maximum error
MaxErr=max(max(abs(un-ue)))
% plot the figure
[X,Y] = pol2cart(The,R);
mesh(X,Y,ue)
set(gca,'fontsize',12)
xlabel('x','fontsize', 16)
ylabel('y','fontsize',16)
zlabel('u','fontsize',16)
view(36,24)

% print -dpng -r600 fdm_polar.png
% print -depsc2 fdm_polar.eps
