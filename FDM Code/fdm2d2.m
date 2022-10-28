% fdm2d2.m
% finite difference method for 2D problem
% -{\Delta}u(x,y)=cos(3*x)*sin(pi*y), (x,y) in (0,pi)x(0,1);
% u(x,0)=u(x,1)=0  in [0,pi]
% u_x(0,y)=u_x(pi,y)=0 in [0,1]
% exact solution: u(x,y)=(9+pi^2)^(-1)*cos(3*x)*sin(pi*y)
clear all; close all;
N=4;   % N=4 8 16 32
% coordinates on the grid
h1=pi/N;  h2=1/N;
x=(0:h1:pi)';
y=(0:h2:1)';
[X,Y]=meshgrid(x,y);
X1=X(2:N,2:N);
Y1=Y(2:N,2:N);
% generate the matrix of RHS
f=cos(3*X1).*sin(pi*Y1);
% constructing the coefficient matrix
e=ones(N-1,1);
C=diag([1/h1^2+2/h2^2, (2/h1^2+2/h2^2)*ones(1,N-3), 1/h1^2+2/h2^2])...
    -1/h1^2*diag(ones(N-2,1),1)-1/h1^2*diag(ones(N-2,1),-1);
D=-1/h2^2*eye(N-1);
% A=kron(eye(N-1),C)+kron(diag(ones(N-2,1),1)+diag(ones(N-2,1),-1),D);
A=kron(eye(N-1),C)+kron(spdiags([e e],[-1 1],N-1,N-1),D);
% solving the linear system
f=f';
un=zeros(N+1,N+1);
un(2:N,2:N)=reshape(A\f(:),N-1,N-1)';
% Neumann boundary condition
un(:,1)=un(:,2);
un(:,end)=un(:,end-1);
ue=1/(9+pi^2)*(cos(3*X)).*(sin(pi*Y));

format long
% value of u and ue in the selected points (i*pi/4,j/4), i,j=1,2,3.
u_select=un(N/4+1:N/4:3*N/4+1,N/4+1:N/4:3*N/4+1)
ue_select=un(N/4+1:N/4:3*N/4+1,N/4+1:N/4:3*N/4+1)

