% LegenSM2_error.m
% Legendre-Galerkin Method for the model equation
% -mu*u''(x)+nu*u'(x)+rho*u(x)=f(x), x in (-1,1),
% boundary condition: u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*xv); 
% RHS: f(x)=mu*kw^2*pi^2*sin(kw*pi*x)+nu*kw*pi*cos(kw*pi*x)+rho*sin(kw*pi*x); 
% Rmk: Use routines lepoly(); legs(); lepolym();
clear, clf
kw=1;
mu=1;
nu=1;
rho=1;
Nvec=[4:2:24];
% Initialization for error
L2_Error=[];  Max_Error=[]; 

for N=Nvec
    [xv,wv]=legs(N+1);       % Legendre-Gauss nodes and weights
    Lm=lepolym(N,xv);        % matrix of Legendre polynomals
    u=sin(kw*pi*xv);         % test function
    f=mu*kw*kw*pi^2*sin(kw*pi*xv)+nu*kw*pi*cos(kw*pi*xv)+rho*sin(kw*pi*xv);  % RHS
    % Calculating coefficients matrix
    S=eye(N-1);              % stiffness matrix
    M=diag(1./(4*[0:N-2]+6))*diag(2./(2*[0:N-2]+1)+2./(2*[0:N-2]+5))...
        -diag(2./(sqrt(4*[0:N-4]+6).*sqrt(4*[0:N-4]+14).*(2*[0:N-4]+5)),2)...
        -diag(2./(sqrt(4*[2:N-2]-2).*sqrt(4*[2:N-2]+6).*(2*[2:N-2]+1)),-2);   % mass matrix
    D=diag(1./(sqrt(2.*[0:N-3]+3).*sqrt(2.*[0:N-3]+5)),1)...
        +diag(-1./(sqrt(2.*[0:N-3]+3).*sqrt(2.*[0:N-3]+5)),-1);     % matrix of derivative term
    A=mu*S+rho*M+nu*D;          % coefficient matrix
    % Solving the linear system
    Pm=diag(1./sqrt(4*[0:N-2]+6))*(Lm(1:end-2,:)-Lm(3:end,:));    % matrix of Phi(x)
    b=Pm*diag(wv)*f; 
    uh=A\b;           % expansion coefficients of u_N(x)
    un=Pm'*uh;        % Coefficiets to points

    L2_error=sqrt(((un-u).^2)'*wv);  % L^2 error
    Max_error=norm(abs(un-u),inf);   % maximum pointwise error 
    L2_Error=[L2_Error;L2_error];
    Max_Error=[Max_Error;Max_error];
end
% Plot L^2 and maximum pointwise error
plot(Nvec,log10(L2_Error),'bo-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(Nvec,log10(Max_Error),'rd-','MarkerFaceColor','w','LineWidth',1)
grid on
legend('L^2 error','L^{\infty} error','location','NorthEast')
% title('Error of Legendre-Galerkin methods','fontsize',12)
set(gca,'fontsize',12)
xlabel('N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)

% sets axis tick and axis limits
xticks(0:5:25)
yticks(-16:2:0)
xlim([0 25])
ylim([-16 0])

% print -dpng -r600  LegenSM2_error.png
