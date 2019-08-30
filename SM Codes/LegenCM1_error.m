% Legendre-collocation method for the model equation:
% -u''(x)+\alpha u(x)=f(x), x in (-1,1);  
% boundary condition: u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x); 
% RHS: f=kw*kw*pi^2*sin(kw*pi*x)+alpha*sin(kw*pi*x); 
% Rmk: Use routines lepoly(); legslb(); legslbdm();
clear all
alpha=1;
kw=10; 
N=32; 
Nvec=[32:2:68]; 
Errv=[];
for N=Nvec
    [x,w]=legslb(N);           % compute LGL nodes and weights
    u=sin(kw*pi*x);            % test solution
    udprime=-kw*kw*pi*pi*sin(kw*pi*x);   
    f=-udprime+alpha*u;    % RHS
    % Setup and solve the collocation system
    D1=legslbdm(N);        % 1st order differentiation matrices
    %D1=legslbdiff(N,x);   % 1st order differentiation matrices
    D2=D1*D1;                 % 2nd order differentiation matrices
    D=(-D2(2:N-1,2:N-1)+alpha*eye(N-2)); % coefficient matrix
    b=f(2:N-1);                 % RHS
    un=D\b;
    un=[0;un;0];               % Solve the system
    
    error=norm(abs(un-u),inf);  % maximum pointwise error
    Errv=[Errv;error];
end;
  plot(Nvec,log10(Errv),'rd-','MarkerFaceColor','w','LineWidth',1.5)
  grid on, xlabel('N','fontsize', 14), ylabel('log10(Error)','fontsize',14)
  title('Convergence of Legendre-collocation method','fontsize',12)
  set(gca,'fontsize',12) 
  
  