% Legendre-Galerkin Method for for the model equation
% -u_xx+u=f in (-1,1) with boundary condition u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x); 
% RHS: f=kw*kw*pi^2*sin(kw*pi*x)+sin(kw*pi*x); 
% Rmk: Use routines lepoly(); legs(); lepolym();
clear all;  clf
kw=10;
Nvec=[32:2:68];    % kw=10
%Nvec=[4:2:22]    % kw=1
Errv=[];                 % Initialization for error
for N=Nvec
    [xv,wv]=legs(N);            % Legendre-Gauss points and weights
    Lm=lepolym(N+1,xv);    % Lm is a Legendre polynomal matrix 
    u=sin(kw*pi*xv);            % test function
    f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv);      % Right-hand-side(RHS)
    % Calculting coefficient matrix
    S=eye(N);         % stiff matrix
    M=diag(1./(4*[0:N-1]+6))*diag(2./(2*[0:N-1]+1)+2./(2*[0:N-1]+5))...
        -diag(2./(sqrt(4*[0:N-3]+6).*sqrt(4*[0:N-3]+14).*(2*[0:N-3]+5)),2)...
        -diag(2./(sqrt(4*[2:N-1]-2).*sqrt(4*[2:N-1]+6).*(2*[2:N-1]+1)),-2);    % mass matrix
    A=S+M;
    % Solving the linear system
    B=diag(1./sqrt(4*[0:N-1]+6))*(Lm(1:end-2,:)-Lm(3:end,:));   % RHS calculative matrix
    b=B*diag(wv)*f;          % Solving RHS
    uh=A\b;                     % expansion coefficients of u_N in terms of the basis
    un=B'*uh;                   % compositing the numerical solution
    
    error=norm(abs(un-u),inf);  % maximum pointwise error 
    Errv=[Errv;error];
end
% Plot the maximum pointwise error
plot(Nvec,log10(Errv),'ro-','MarkerFaceColor','w','LineWidth',1.5)
grid on, 
xlabel('N','fontsize', 14), ylabel('log10(Error)','fontsize',14)
title('Round-off error of Legendre-Galerkin methods','fontsize',12)
set(gca,'fontsize',12)


print -dpng -r600  LegenSM1_error.png
