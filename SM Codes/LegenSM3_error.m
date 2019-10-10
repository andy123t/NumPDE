% Legendre-Spectral Method for 1D elliptic problem
% -u_yy+u=f in [0,1] with boundary condition: u(0)=1, u'(1)=0;
% exact solution: u=(1-y)^2*exp(y);   RHS: f=(2-4*y)*exp(y);
% Converted : -4U_xx+U=F in [-1,1]
% boundary condition: U(-1)=0, U'(1)=0;
% exact solution: U=(1/2-1/2*x)^2*exp(1/2*x+1/2)-1;
% RHS:  F=-2*x*exp(1/2*x+1/2)-1.
clear all;  clf
Nvec=2:16;  
Errv=[];  condnv=[];                % Initialization for error and condition number
for N=Nvec
    [xv,wv]=legs(N);               % xv and wv are Legendre-Gauss points and weights
    Lm=lepolym(N+1,xv);      % Lm is a Legendre polynomal matrix 
    yv=1/2*(xv+1);                 % variable substitution 
    U=(1-yv).^2.*exp(yv)-1;    % test function
    F=(2-4*yv).*exp(yv)-1;       % RHS in [0,1] 
    
    % Calculting coefficient matrix
    e1=0:N-1;  e2=0:N-2;  e3=0:N-3;
    S=diag( (4*e1+6).*(e1+1).^2./(e1+2).^2 );     % stiff matrix
    M=diag(2./(2*e1+1)+2*(2*e1+3)./(e1+2).^4+2*((e1+1)./(e1+2)).^4./(2*e1+5))...
        +diag( 2./(e2+2).^2-2*(e2+1).^2./((e2+2).^2.*(e2+3).^2) , 1 )...
        +diag(2./(e2+2).^2-2*(e2+1).^2./((e2+2).^2.*(e2+3).^2),-1)...
        +diag( -2*(e3+1).^2./((2*e3+5).*(e3+2).^2) , 2 )...
        +diag(-2*(e3+1).^2./((2*e3+5).*(e3+2).^2),-2);  % mass matrix
    A=4*S+M;
    
    % Solving the linear system
    B=(Lm(1:end-2,:)+diag((2*e1+3)./(e1+2).^2)*Lm(2:end-1,:)-diag((e1+1).^2./(e1+2).^2)*Lm(3:end,:));
    b=B*diag(wv)*F;          % Solving RHS
    Uh=A\b;                      % expansion coefficients of u_N in terms of the basis
    Un=B'*Uh;                   % compositing the numerical solution
    
    error=norm(abs(Un-U),2);  % L^2 error 
    Errv=[Errv;error];
    condnv=[condnv,cond(A)];   % condition number of A
end
% Plot the maximum pointwise error
plot(Nvec,log10(Errv),'s-','color',[0 0.5 0],'MarkerFaceColor','w','LineWidth',1.5)
grid on, 
xlabel('N','fontsize', 14), ylabel('log_{10}(Error)','fontsize',14)
title('L^2 error of Legendre-Galerkin method','fontsize',12)
set(gca,'fontsize',12)


print -dpng -r600  LegenSM3_error.png