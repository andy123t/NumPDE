% LegenSM1_error.m
% Legendre-Galerkin Method for the model equation
% -u_xx+u=f in (-1,1) with boundary condition u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x); 
% RHS: f=kw*kw*pi^2*sin(kw*pi*x)+sin(kw*pi*x); 
% Rmk: Use routines lepoly(); legs(); lepolym();
clear all;  clf
kw=10;
Nvec=[32:2:76];

L2_Error=[]; Max_Error=[];    % initialization error

for N=Nvec
    [xv,wv]=legs(N+1);        % Legendre-Gauss nodes and weights
    Lm=lepolym(N,xv);         % matrix of Legendre polynomals
    u=sin(kw*pi*xv);          % test function
    f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv);   % Right-hand-side(RHS)
    % Calculting coefficient matrix
    S=eye(N-1);               % stiff matrix
    M=diag(1./(4*[0:N-2]+6))*diag(2./(2*[0:N-2]+1)+2./(2*[0:N-2]+5))...
        -diag(2./(sqrt(4*[0:N-4]+6).*sqrt(4*[0:N-4]+14).*(2*[0:N-4]+5)),2)...
        -diag(2./(sqrt(4*[2:N-2]-2).*sqrt(4*[2:N-2]+6).*(2*[2:N-2]+1)),-2);    % mass matrix
    A=S+M;
    % Solving the linear system
    Pm=diag(1./sqrt(4*[0:N-2]+6))*(Lm(1:end-2,:)-Lm(3:end,:));   % matrix of Phi(x)
    b=Pm*diag(wv)*f;          % solving RHS
    uh=A\b;                   % expansion coefficients of u_N(x)
    un=Pm'*uh;                % compositing the numerical solution
    
    L2_error=norm(abs(un-u),2);     % L^2 error
    Max_error=norm(abs(un-u),inf);  % maximum pointwise error 
    L2_Error=[L2_Error;L2_error];
    Max_Error=[Max_Error;Max_error];
end
% Plot L^2 and maximum pointwise error
plot(Nvec,log10(L2_Error),'ro-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(Nvec,log10(Max_Error),'md-','MarkerFaceColor','w','LineWidth',1)
grid on
legend('L^2 error','L^{\infty} error','location','NorthEast')
% title('Error of Legendre-Galerkin methods','fontsize',12)
set(gca,'fontsize',12)
xlabel('N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)

% sets axis tick and axis limits
xticks(30:10:80)
yticks(-15:3:0)
xlim([30 80])
ylim([-15 0])

% print -dpng -r600  LegenSM1_error.png
