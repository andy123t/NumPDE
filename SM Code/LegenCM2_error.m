% Legendre-collocation Method for the model equation: 
% -u''(x)+u'(x)+u(x)=f(x),  x in (-1,1);
% % boundary condition: u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x); 
% RHS: f(x)=kw^2*pi^2*sin(kw*pi*x)+sin(kw*pi*x)+kw*pi*cos(kw*pi*x); 
% Rmk: Use routines lepoly(); legslb(); legslbdm(); 
clear all;  clf
kw=10;
Nv=[32:2:68]; 
Errv=[];
for N=Nv
    [xv,wv]=legslb(N);          % compute LGL nodes and weights
    u=sin(kw*pi*xv);            % test function
    f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv)+kw*pi*cos(kw*pi*xv);  % RHS
    % Setup and solve the collocation system
    D1=legslbdm(N);   %1st order differentiation matrix
    D2=D1*D1;            % 2nd order differentiation matrix
    D=-D2(2:N-1,2:N-1)+D1(2:N-1,2:N-1)+eye(N-2);    % coefficient matrix
    b=f(2:N-1);   % RHS
    un=D\b;
    un=[0;un;0];  % Solve the system 
    
    error=norm(abs(un-u),inf); 
    Errv=[Errv;error];
end
% Plot the maximum pointwise error 
plot(Nv,log10(Errv),'md-','MarkerFaceColor','w','LineWidth',1.5)
grid on, 
xlabel('N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)
title('L^{\infty} error of Legendre-collocation method','fontsize',12)


print -dpng -r600  LegenCM2_error.png
