% ChebSM1_error.m
% Chebyshev-Galerkin Method for the model equation
% -u_xx+u=f in (-1,1) with boundary condition u(-1)=u(1)=0;
% exact solution: u=sin(kw*pi*x);
% RHS: f=kw*kw*pi^2*sin(kw*pi*x)+sin(kw*pi*x);
% Rmk: Use routines japoly(); jags(); japolym();
clear, clf
kw=10;
Nvec=32:2:76;
% Initialization for error
L2_Err=[];  Max_Err=[];

for N=Nvec
    [xv,wv]=jags(N+1,-1/2,-1/2);               % Chebyshev-Gauss points and weights
    Cm=japolym(N,-1/2,-1/2,xv)./japolym(N,-1/2,-1/2,1);  % matrix of Chebyshev polynomal
    u=sin(kw*pi*xv);                           % test function
    f=kw*kw*pi^2*sin(kw*pi*xv)+sin(kw*pi*xv);  % RHS of test function
    % Calculting coefficient matrix
    S=zeros(N-1);    % stiff matrix
    for k=1:N-1
        for j=1:N-1
            if k==j
                S(k,j)=2*pi*k*(k+1);
            elseif (k<j && mod(j-k,2)==0)
                S(k,j)=4*pi*k;
            else
                S(k,j)=0;
            end
        end
    end
    
    M=diag([3/2*pi, pi*ones(1,N-2)])+diag(-pi/2*ones(1,N-3),2)...
        +diag(-pi/2*ones(1,N-3),-2);    % mass matrix
    A=S+M;
    
    % Solving the linear system
    Pm=Cm(1:end-2,:)-Cm(3:end,:);   % matrix of Phi(x)
    b=Pm*diag(wv)*f;          % Solving RHS
    uh=A\b;                   % expansion coefficients of u_N(x)
    un=Pm'*uh;                % compositing the numerical solution
    
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
% title('L^2 error of Chebyshev-Galerkin methods','fontsize',12)
set(gca,'fontsize',12)
xlabel('N','fontsize', 14), ylabel('log_{10}Error','fontsize',14)

% sets axis tick and axis limits
xticks(30:10:80)
yticks(-14:2:0)
xlim([30 80])
ylim([-14 0])

% print -dpng -r600  ChebSM1_error.png
