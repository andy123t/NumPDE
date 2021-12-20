% FEM1D.m
% finite element method for 1D elliptic problem
% -u_xx+u=f in (0,1)
% boundary condition: u(0)=u(1)=0;
% exact solution: u=x*(1-x)*sin(x)
% RHS: f=(4*x-2).*cos(x)+(2+2*x-2*x^2).*sin(x);
clear all; close all;
Num=[16 32 64 128 256 512];  % Number of partitions
Err=[];  DOF=[];
for j=1:length(Num)
    N=Num(j);  h=1/N;  x=0:h:1;
    % The global node number corresponds to the element local node number
    M=[1:N;2:N+1];
    [xv,wv]=jags(2,0,0); % nodes and weights of gauss quadrature
    
    K=zeros(N+1);        % global stiffness matrix
    F=zeros(N+1,1);      % RHS load vector
    for i=1:N   % loop for each element
        K(M(1,i),M(1,i))=K(M(1,i),M(1,i))+((h/2)*(((1/4)*(2/h)^2+((1-xv)/2).^2)))'*wv;
        K(M(1,i),M(2,i))=K(M(1,i),M(2,i))+((h/2)*((-1/4)*(2/h)^2+((1-xv)/2).*((1+xv)/2)))'*wv;
        K(M(2,i),M(1,i))=K(M(2,i),M(1,i))+((h/2)*((-1/4)*(2/h)^2+((1-xv)/2).*((1+xv)/2)))'*wv;
        K(M(2,i),M(2,i))=K(M(2,i),M(2,i))+((h/2)*(((1/4)*(2/h)^2+((1+xv)/2).^2)))'*wv;
        
        t=h*xv/2+(x(i+1)+x(i))/2;
        F(M(1,i))=F(M(1,i))+(h/2*((1-xv)/2).*((4*t-2).*cos(t)+(2+2*t-2*t.^2).*sin(t)))'*wv;
        F(M(2,i))=F(M(2,i))+(h/2*((1+xv)/2).*((4*t-2).*cos(t)+(2+2*t-2*t.^2).*sin(t)))'*wv;
    end
    % Dirichlet boundary condition
    K(1,: )=zeros(1,N+1);
    K(:,1)=zeros(1,N+1);
    K(N+1,: )=zeros(1,N+1);
    K(:, N+1)=zeros(1,N+1);
    K(1,1)=1;  K(N+1,N+1)=1;
    F(1)=0;    F(N+1)=0;

    U=K\F;     % numerical solution at the value of the node
    err=max(abs(U'-x.*(1-x).*sin(x)));  % node error
    doff=N+1;  % degrees of freedom, number of unknowns
    Err=[Err, err];
    DOF=[DOF, doff];
end
plot(log10(DOF),log10(Err),'ro-','MarkerFaceColor','w','LineWidth',1),
hold on
plot(log10(DOF),log10(DOF.^(-2)),'--')
grid on
%title('Convergence of Finite Element Method','fontsize',14)
set(gca,'fontsize',14)
xlabel('log_{10}N','fontsize',14), ylabel('log_{10}Error','fontsize',14),

% print -dpng -r600 FEM1D.png
% print -depsc2 FEM1D.png
