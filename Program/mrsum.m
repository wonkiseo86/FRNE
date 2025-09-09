function rs = mrsum(x,y)

% Matrix Riemann Sum(MRSUM)
% input x : n*1 vector
%       y : n*m matrix
%
% output rs : m*1 vector

z  = [x y] ;
n  = length(z) ;
nc = size(y,2) ;

t  = sortrows(z) ;
x0 = t(:,1);
y0 = t(:,2:nc+1) ;
dx = x0(2:n,:)-x0(1:n-1,:) ;
y1 = y0(1:n-1,:) ;
y2 = y0(2:n, :) ;

s1 = y1'*dx ;
s2 = y2'*dx ;

rs = (s1+s2)/2 ;
end