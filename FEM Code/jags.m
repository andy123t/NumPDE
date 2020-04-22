% x=jags(n,alp,bet) returns n Jacobi-Gauss points with parameter (alp,bet)
%                   by using eigenvalue method for nodes.
% [x,w]=jags(n,alp,bet) also returns the weights (stored in w). 
% See Page 84 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
% Last modified on September 4, 2011

function [varargout]=jags(n,alp,bet)

if n<=0, disp('Input n >=1'); varargout{1}='Wrong input';  return; end


apb=alp+bet;
if n==1, 
    varargout{1}=(bet-alp)/(apb+2);
    varargout{2}=exp((apb+1)*log(2)+gammaln(alp+1)+gammaln(bet+1)-gammaln(apb+1)); 
    return; 
end;


if alp==bet, 
  dm=zeros(1,n);
else 
 dm=[(bet-alp)/(2+apb), (bet^2-alp^2)./((2*[1:n-1]+apb).*(2*[1:n-1]+apb+2))];            %  Create maindiagnal
end;


id=[2:n-1];                                              % Indices     
ds=[4*(1+alp)*(1+bet)/((2+apb)^2*(3+apb)),4*id.*(id+alp).*(id+bet).*(id+apb)./((2*id+apb-1).*(2*id+apb).^2.*(2*id+apb+1))];   %  Create subdiagonals
ds=sqrt(ds);

J =diag(dm)+diag(ds,1)+diag(ds,-1);      %  Create Jacobi matrix
x= sort(eig(sparse(J)));                 %  Compute eigenvalues

varargout{1}=x;                        % Return nodes 

if nargout==1, return; end;

 gn=(apb+1)*log(2)+gammaln(n+alp+1)+gammaln(n+bet+1)-gammaln(n+1)-gammaln(n+apb+1);
 gn=exp(gn);                         % Constant in the weight expression (refer to (3.132b)of Shen Tang and Wang's book)
 [dy,y]=japoly(n,alp,bet,x);         % Compute derivative of Jacobi polynomial of degree n 
                                     % at the nodes 
 varargout{2}=gn./((1-x.^2).*dy.^2);  % Compute weights by using weight expression. 

return;


 
   
   
 
  
  
