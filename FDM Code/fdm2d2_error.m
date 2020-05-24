% fdm2d2_error.m 
% finite difference method for 2D problem
% -\Delta u=cos(3*x)*sin(pi*y)  in (0,pi)x(0,1)
% u(x,0)=u(x,1)=0  in [0,pi]
% u_x(0,y)=u_x(pi,y)=0 in [0,1]
% exact solution: ue=(9+pi^2)^(-1)*cos(3*x)*sin(pi*y)
clear all; close all; clf
Nvec=2.^[2:7];
Error=[];
for N=Nvec
    % coordinates on the grid
    h1=pi/N;  h2=1/N;
    x=[0:h1:pi]'; 
    y=[0:h2:1]';
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
    %A=kron(eye(N-1),C)+kron(diag(ones(N-2,1),1)+diag(ones(N-2,1),-1),D);
    A=kron(eye(N-1),C)+kron(spdiags([e e],[-1 1],N-1,N-1),D);
    % solving the linear system
    f=f'; 
    u=zeros(N+1,N+1);
    u(2:N,2:N)=reshape(A\f(:),N-1,N-1)';
    % Neumann boundary condition
    u(:,1)=u(:,2); 
    u(:,end)=u(:,end-1);
    ue=1/(9+pi^2)*(cos(3*X)).*(sin(pi*Y));
    error=max(max(abs(u-ue)));   %maximum error
    Error=[Error,error];
end
plot(log10(Nvec),log10(Error),'ro-','MarkerFaceColor','w','LineWidth',1)
hold on,
plot(log10(Nvec), log10(Nvec.^(-1)), '--')
grid on,
%title('Convergence of Finite Difference Method','fontsize',14)
set(gca,'fontsize',14)
xlabel('log_{10}N','fontsize', 14), ylabel('log_{10}Error','fontsize',14),

% add annotation of slope
ax = [0.64 0.60];
ay = [0.69 0.64];
annotation('textarrow',ax,ay,'String','slope = -1 ','fontsize',14)

% computating convergence order
for i=1:length(Nvec)-1
    order(i)=-log(Error(i)/Error(i+1))/(log(Nvec(i)/Nvec(i+1)));
end
Error
order

% print -dpng -r600  fdm2d2_error.png
