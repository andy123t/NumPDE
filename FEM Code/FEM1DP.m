% FEM1DP.m
% finite element method for 1D elliptic problem
% -u_xx+u=f in (0,1)
% boundary condition: u(0)=u(1)=0;
% exact solution: u=x*(1-x)*sin(x);
% right hand function: f=(4*x-2).*cos(x)+(2+2*x-2*x^2).*sin(x)
clear all; close all;
Num=[16 32 64 128 256 512];  % Number of partitions
node_Err=[];  L2_Err=[];  H1_Err=[];  DOF=[];
for k=1:length(Num)
    N=Num(k);  h=1/N;  x=0:h:1;
    
    M=[1:N;2:N+1];        % information matrix
    % the global node indices of the mesh nodes of all the mesh elements
    [xv,wv]=jags(3,0,0);  % nodes and weights of Gauss quadrature

    K=zeros(N+1);         % global stiffness matrix
    F=zeros(N+1,1);       % RHS load vector

    for i=1:N             % loop for each element
        K(M(1,i),M(1,i))=K(M(1,i),M(1,i))+((h/2)*(((1/4)*(2/h)^2+((1-xv)/2).^2)))'*wv;
        K(M(1,i),M(2,i))=K(M(1,i),M(2,i))+((h/2)*((-1/4)*(2/h)^2+((1-xv)/2).*((1+xv)/2)))'*wv;
        K(M(2,i),M(1,i))=K(M(2,i),M(1,i))+((h/2)*((-1/4)*(2/h)^2+((1-xv)/2).*((1+xv)/2)))'*wv;
        K(M(2,i),M(2,i))=K(M(2,i),M(2,i))+((h/2)*(((1/4)*(2/h)^2+((1+xv)/2).^2)))'*wv;
        
        t=h*xv/2+(x(i+1)+x(i))/2;
        F(M(1,i))=F(M(1,i))+(h/2*((1-xv)/2).*((4*t-2).*cos(t)+(2+2*t-2*t.^2).*sin(t)))'*wv;
        F(M(2,i))=F(M(2,i))+(h/2*((1+xv)/2).*((4*t-2).*cos(t)+(2+2*t-2*t.^2).*sin(t)))'*wv;
    end
    
    % Handling Dirichlet boundary condition
    K(1,: )=zeros(1,N+1);
    K(:,1)=zeros(1,N+1);
    K(N+1,: )=zeros(1,N+1);
    K(:, N+1)=zeros(1,N+1);
    K(1,1)=1;  K(N+1,N+1)=1;
    F(1)=0;    F(N+1)=0;

    U=K\F;       % numerical solution at the value of the nodes
    node_error=max(abs(U'-x.*(1-x).*sin(x)));  % node error
    for i=1:N
        tt=h*xv/2+(x(i+1)+x(i))/2;
        % value of finite element solution at Gauss nodes
        uh=U(i)*(1-xv)/2+U(i+1)*(1+xv)/2;
        % derivative value of finite element solution at Gauss nodes
        duh=-U(i)/2+U(i+1)/2;
        L2_err(i)=h/2*((tt.*(1-tt).*sin(tt)-uh).^2)'*wv;
        % the square of the L2 error of the i-th interval
        H1_err(i)=h/2*((sin(tt)-2*tt.*sin(tt)+tt.*(1-tt).*cos(tt)-duh*2/h).^2)'*wv;
        % the square of the H1 semi-norm error of the i-th interval
    end
    node_Err=[node_Err, node_error];
    L2_Err=[L2_Err, sqrt(sum(L2_err))];
    H1_Err=[H1_Err, sqrt(sum(L2_err)+sum(H1_err))];
    doff=N+1;    % degrees of freedom, number of unknowns
    DOF=[DOF, doff];
end
loglog(DOF,node_Err,'r+-','LineWidth',1)
hold on
loglog(DOF,L2_Err,'bo-','MarkerFaceColor','w','LineWidth',1)
hold on
loglog(DOF,H1_Err,'m*-','LineWidth',1)
hold on
grid on
legend('L^2 error','L^{\infty} error','H^1 error','location','SouthWest')
% title('Convergence of finite element method','fontsize',14)
set(gca,'FontName','Times New Roman','FontSize',12)
xlabel('Number of degrees of freedom','fontsize', 14),
ylabel('Error','fontsize',14)

% sets axis tick and axis limits
xticks(10.^(1:3))
yticks(10.^(-8:-1))
xlim([10 10^3])
ylim([10^(-8) 10^(-1)])

% calculating of convergence order
for k=1:length(Num)-1
    node_order(k)=log(node_Err(k)/node_Err(k+1))/(log(DOF(k)/DOF(k+1)));
    L2_order(k)=log(L2_Err(k)/L2_Err(k+1))/(log(DOF(k)/DOF(k+1)));
    H1_order(k)=log(H1_Err(k)/H1_Err(k+1))/(log(DOF(k)/DOF(k+1)));
end
node_order
L2_order
H1_order

% print -dpng -r600 FEM1DP.png
% print -depsc2 FEM1DP.eps
