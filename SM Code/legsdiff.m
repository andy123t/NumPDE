% D=legsdiff(n,x) returns the first-order differentiation matrix of size
% n by n, associated with the Legendre-Gauss points x, which may be computed by 
% x=legs(n) or x=legsndm(n). 
% See Page 110 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepoly() 
% Last modified on August 31, 2011

function D=legsdiff(n,x)
if n==0, D=[]; return; end;
xx=x;[dy,y]=lepoly(n,xx); nx=size(x); 
if nx(2)>nx(1), dy=dy'; y=y'; xx=x'; end;  %% dy is a column vector of L_{n}'(x_k)
  D=(xx./dy)*dy'-(1./dy)*(xx.*dy)';  %% compute L_{n}'(x_j) (x_k-x_j)/L_{n}'(x_k);     
                                 % 1/d_{kj} for k not= j (see (3.203)) 
  D=D+eye(n);                    % add the identity matrix so that 1./D can be operated                                     
  D=1./D; 
  D=D-eye(n); D=D+diag(xx./(1-xx.^2));  % update the diagonal entries  
  return; 
 