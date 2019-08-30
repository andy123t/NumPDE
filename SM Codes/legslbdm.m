 function d=legslbdm(n);

%  The function d=legslbdm(n) computes the first-order Differentiation Matrix (N-by-N) w.r.t 
%  Legndre-Gauss-Lobatto 
%  Rmk: Use routines lepoly() and legslb();
%  Last Modified: 02/06/2005.


if n<2, disp('n>2 is necessary'); return; end;

x=legslb(n);   %Compute the nodes 
y=lepoly(n-1,x);
%% Evaluate the derivative matrix 
dd=(x./y)*y'-(1./y)*(x.*y)';
dd=dd+eye(size(dd));d=1./dd;
d=d-eye(size(d));d(1,1)=n*(n-1)/4; d(n,n)=-d(1,1);

% Test the case 
 %y=sin(x)+1;dy=cos(x);  d*y-dy
return; 



