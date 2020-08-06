%Y=japolym(n,alp,bet,x) computes the Jacobi polynomials of degree up to n with parameter (alp,bet) at a
%   vector-valued x
% [DY,Y]=japolym(n,alp,bet, x) also returns the values of the 1st-order 
%   derivatives stored in DY
% Note: Y (and likewise for DY) saves J_0^{alp,bet}(x), J_1^{alp,bet}(x), ...., J_n^{alp,bet}_n(x) by rows
% i.e., J_k^{alp,bet}(x) is the (k+1)th row of the matrix Y (or DY)
% See Page 74 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
% 
% Last modified on September 2, 2011    


function [varargout]=japolym(n,alp,bet,x)
  apb=alp+bet;     

dim=size(x); xx=x; if dim(1)>dim(2), xx=xx'; end; % xx is a row-vector
if nargout==1,
     if n==0, varargout{1}=ones(size(xx));  return; end;
     if n==1, varargout{1}=[ones(size(xx)); 0.5*(alp-bet+(apb+2)*xx)]; return; end;

     polylst=ones(size(xx));	  poly=0.5*(alp-bet+(apb+2)*xx);
     Y=[polylst;poly];
     
   for k=2:n,
	  a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	  a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	  b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	  a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	  polyn=((a2+a3*xx).*poly-a4*polylst)/a1;  % See (3.110) and (3.111) 
      polylst=poly; poly=polyn;	Y=[Y;polyn];
   end;
      varargout{1}=Y; return;
end;

   
if n==0, varargout{2}=ones(size(xx)); varargout{1}=zeros(size(xx)); return; end;
if n==1, varargout{2}=[ones(size(xx)); 0.5*(alp-bet+(apb+2)*xx)];
    varargout{1}=[zeros(size(xx));0.5*(apb+2)*ones(size(xx))]; return; end;

   polylst=ones(size(xx));         pderlst=zeros(size(xx));
   poly=0.5*(alp-bet+(apb+2.)*xx); pder=0.5*(apb+2.)*ones(size(xx));
   Y=[polylst;poly]; DY=[pderlst;pder];
   for k=2:n,
	  a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	  a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	  b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	  a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	  polyn=((a2+a3*xx).*poly-a4*polylst)/a1;  % See (3.110) and (3.111) 
	  pdern=((a2+a3*xx).*pder-a4*pderlst+a3*poly)/a1;	  	  
	  polylst=poly; poly=polyn; Y=[Y;polyn];
	  pderlst=pder; pder=pdern; DY=[DY;pdern];
   end;
      varargout{2}=Y;
      varargout{1}=DY;
   return;


    

  
     
	
