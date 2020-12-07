% LegenCM2_error.m
% Legendre-collocation Method for the model equation: 
% -mu*u''(x)+nu*u'(x)+rho*u(x)=f(x),  x in (-1,1);
% % boundary condition: u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x); 
% RHS: f(x)=mu*kw^2*pi^2*sin(kw*pi*x)+nu*kw*pi*cos(kw*pi*x)+rho*sin(kw*pi*x); 
% Rmk: Use routines lepoly(); legslb(); legslbdm(); 
clear, clf
kw=10;
nu=1;
mu=1;
rho=1;
Nvec=32:2:76;
% Initialization for error
L2_Err=[];  Max_Err=[];
for N=Nvec
    [xv,wv]=legslb(N);          % compute LGL nodes and weights
    u=sin(kw*pi*xv);            % test function
    f=mu*kw*kw*pi^2*sin(kw*pi*xv)+nu*kw*pi*cos(kw*pi*xv)+rho*sin(kw*pi*xv);  % RHS
    % Setup and solve the collocation system
    D1=legslbdiff(N,xv);  % 1st order differentiation matrix
    %D1=legslbdm(N);      % 1st order differentiation matrix
    D2=D1*D1;             % 2nd order differentiation matrix
    % Compositing the coefficient matrix
    D=-mu*D2(2:N-1,2:N-1)+nu*D1(2:N-1,2:N-1)+rho*eye(N-2);
    b=f(2:N-1);   % RHS
    un=D\b;
    un=[0;un;0];  % Solve the system 
    
    L2_error=sqrt(((un-u).^2)'*wv);  % L^2 error
    Max_error=norm(abs(un-u),inf);   % maximum pointwise error 
    L2_Err=[L2_Err;L2_error];
    Max_Err=[Max_Err;Max_error];
end
% Plot L^2 and maximum pointwise error
plot(Nvec,log10(L2_Err),'bo-','MarkerFaceColor','w','LineWidth',1)
hold on
plot(Nvec,log10(Max_Err),'rd-','MarkerFaceColor','w','LineWidth',1)
grid on
legend('L^2 error','L^{\infty} error','location','NorthEast')
set(gca,'fontsize',12)
xlabel('N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)
%title('Convergence of Legendre-collocation method','fontsize',12)

% sets axis tick and axis limits
xticks(30:10:80)
yticks(-15:3:0)
xlim([30 80])
ylim([-15 0])

% print -dpng -r600  LegenCM2_error.png
