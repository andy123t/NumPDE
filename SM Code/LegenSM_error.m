% LegenSM_error.m
% Legendre-Galerkin Method for the model equation
% -u_xx+u=f in (-1,1) with boundary condition u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x);
% RHS: f=kw*kw*pi^2*sin(kw*pi*x)+sin(kw*pi*x);
% Rmk: Use routines lepoly(); legs(); lepolym();
clear all;  clf
kw=1;
Nvec=[4:2:28];              % Degree of polynomials
L2_Error=[]; Max_Error=[];  % Initialization for error
for N=Nvec
    [xv,wv]=legs(N+1);         % Legendre-Gauss nodes and weights
    Lm=lepolym(N,xv);          % matrix of Legendre polynomals
    u=sin(kw*pi*xv);           % test function
    f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv);  % Right-hand-side(RHS)
    % Calculting coefficient matrix
    S=diag(4*[0:N-2]+6);       % stiff matrix
    M=diag(2./(2*[0:N-2]+1)+2./(2*[0:N-2]+5))...
        -diag(2./(2*[0:N-4]+5),2)...
        -diag(2./(2*[2:N-2]+1),-2);  % mass matrix
    A=S+M;
    % Solving the linear system
    Pm=Lm(1:end-2,:)-Lm(3:end,:);  % matrix of Phi(x)
    b=Pm*diag(wv)*f;               % Solving RHS
    uh=A\b;                        % expansion coefficients of u_N(x)
    un=Pm'*uh;                     % compositing the numerical solution
    
    Max_error=max(abs(un-u));      % maximum pointwise error
    L2_error=norm(abs(un-u),2);    % L^2 error
    Max_Error=[Max_Error;Max_error];
    L2_Error=[L2_Error;L2_error];
end
% Plot the L^2 and maximum pointwise error
plot(Nvec,log10(L2_Error),'bo-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(Nvec,log10(Max_Error),'rd-','MarkerFaceColor','w','LineWidth',1)
grid on,
legend('L^2 error','L^{\infty} error')
% title('Error of Legendre-Galerkin methods','fontsize',12)
set(gca,'fontsize',12)
xlabel('N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)

% print -depsc2 LegenSM_error.eps
% print -dpng -r600  LegenSM_error.png
