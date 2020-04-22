%y=japoly(n,alp,bet,x) computes the Jacobi polynomial of degree n with parameter (alp,bet) at a
%   vector-valued x
% [dy,y]=japoly(n,alp,bet, x) also returns the values of the 1st-order 
%   derivative stored in dy
% See Page 74 of the book: J. Shen, T. Tang and L. Wang, Spectral Methods:
%  Algorithms, Analysis and Applications, Springer Series in Compuational
%  Mathematics, 41, Springer, 2011. 
% 
% Last modified on September 2, 2011    


function [varargout]=japoly(n,alp,bet,x)
  apb=alp+bet;     
   
if nargout==1,
     if n==0, varargout{1}=ones(size(x));  return; end;
     if n==1, varargout{1}=0.5*(alp-bet+(apb+2)*x); return; end;

     polylst=ones(size(x));	
     poly=0.5*(alp-bet+(apb+2)*x);
     
   for k=2:n,
	  a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	  a2=(2.0*k+apb-1.)*(alp^2-bet^2);  %(b3+1)*(alp^2-bet^2);--Jiao remark
	  b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	  a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);  %2.0*(k+alp-1.)*(k+bet-1.)*(b3+2)--Jiao remark
	  polyn=((a2+a3*x).*poly-a4*polylst)/a1;  % See (3.110) and (3.111) 
      polylst=poly; poly=polyn;	
   end;
      varargout{1}=polyn; return;
end;

   
if n==0, varargout{2}=ones(size(x)); varargout{1}=zeros(size(x)); return; end;
if n==1, varargout{2}=0.5*(alp-bet+(apb+2)*x); varargout{1}=0.5*(apb+2)*ones(size(x)); return; end;

   polylst=ones(size(x));         pderlst=zeros(size(x));
   poly=0.5*(alp-bet+(apb+2.)*x); pder=0.5*(apb+2.)*ones(size(x));
   
   for k=2:n,
	  a1=2.0*k*(k+apb)*(2.0*k+apb-2.);
	  a2=(2.0*k+apb-1.)*(alp^2-bet^2);
	  b3=(2.0*k+apb-2.); a3=b3*(b3+1.)*(b3+2.);
	  a4=2.0*(k+alp-1.)*(k+bet-1.)*(2.0*k+apb);
	  polyn=((a2+a3*x).*poly-a4*polylst)/a1;  % See (3.110) and (3.111) 
	  pdern=((a2+a3*x).*pder-a4*pderlst+a3*poly)/a1;	  	  
	  polylst=poly; poly=polyn;
	  pderlst=pder; pder=pdern;
   end;
      varargout{2}=polyn;
      varargout{1}=pdern;
   return;


    

  
     
	
