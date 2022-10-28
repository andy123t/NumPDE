 function [varargout]=legs(n)

%  x=legs(n) returns n Legendre-Gauss points arranged in ascending order
%  [x,w]= legs(n) returns n Legendre-Gauss points and weights
%  Newton iteration method is used for computing nodes
%  See Page 99 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
%  Use the function: lepoly() 
%  Last modified on August 30, 2011

 % Compute the initial guess of the interior LGL points
  thetak=(4*[1:n]-1)*pi/(4*n+2);
  ze=-(1-(n-1)/(8*n^3)-(39-28./sin(thetak).^2)/(384*n^4)).*cos(thetak);
  ep=eps*10;                            % error tolerance for stopping iteration
  ze1=ze+ep+1;
 
 while max(abs(ze1-ze))>=ep,            % Newton's iteration procedure
      ze1=ze;
      [dy,y]=lepoly(n,ze);
      ze=ze-y./dy;  % see Page 99 of the book
 end;                                   % around 6 iterations are required for n=100
    varargout{1}=ze';
 if nargout==1, return; end;
   
 % Use the weight expression (3.178) to compute the weights
  varargout{2}=(2./((1-ze.^2).*dy.^2))';
  




