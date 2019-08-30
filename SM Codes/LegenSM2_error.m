% Legendre-Galerkin Method for the model equation
% -u''(x)+u'(x)+u(x)=f(x), x in (-1,1),
% boundary condition: u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*xv); 
% RHS: f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv); 
% Rmk: Use routines lepoly(); legs(); lepolym();
clear all
kw=10;
Nvec=[32:2:68];
Errv=[];
for N=Nvec
    [xv,wv]=legs(N);             % Legendre-Gauss points and weights
    Lm=lepolym(N+1,xv);    % Lm is a Legendre polynomal matrix 
    u=sin(kw*pi*xv);            % test function
    f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv)+kw*pi*cos(kw*pi*xv);  % RHS
    % Calculating coefficients matrix
    S=eye(N);   % stiffness matrix
    M=diag(1./(4*[0:N-1]+6))*diag(2./(2*[0:N-1]+1)+2./(2*[0:N-1]+5))...
        -diag(2./(sqrt(4*[0:N-3]+6).*sqrt(4*[0:N-3]+14).*(2*[0:N-3]+5)),2)...
        -diag(2./(sqrt(4*[2:N-1]-2).*sqrt(4*[2:N-1]+6).*(2*[2:N-1]+1)),-2);  % mass matrix
    D=diag(1./(sqrt(2.*[0:N-2]+3).*sqrt(2.*[0:N-2]+5)),1)...
        +diag(-1./(sqrt(2.*[0:N-2]+3).*sqrt(2.*[0:N-2]+5)),-1);          % matrix derived from u'(x)
    A=S+M+D;       % Coefficient matrix
    % Solving the linear system
    B=diag(1./sqrt(4*[0:N-1]+6))*(Lm(1:end-2,:)-Lm(3:end,:));   % RHS calculative matrix
    b=B*diag(wv)*f; 
    uh=A\b;          % expansion coefficients of u_N in terms of the basis
    un=B'*uh;        % Coefficiets to points
    error=norm(abs(un-u),inf);   % maximum pointwise error
    Errv=[Errv;error]; 
end
% Plot the maximum pointwise error
plot(Nvec,log10(Errv),'mo-','MarkerFaceColor','w','LineWidth',1.5)
grid on, xlabel('N','fontsize', 14), ylabel('log10(Error)','fontsize',14) 
title('Round-off error of Legendre-Galerkin methods','fontsize',12)
set(gca,'fontsize',12)

