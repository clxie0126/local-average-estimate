function b = xia(data,h)
%   Xia et al. (2004) for y=x'a(u)+z'b+e
%   semi-local least squares estimator
%   data=[u x z y]
%   local linear fit for coefficient functions
%  

q=1; %the dimension of z

[n, m]=size(data);

u=data(:,1);
x=data(:,2:m-1-q);
z=data(:,m-q:m-1);
Y=data(:,m);

p=m-q-2;
A=kron(eye(n),[x,diag(u./h)*x])-kron(diag(u./h),[zeros(size(x)),x]);
omega=[A,kron(ones(n,1),z)];
du=kron(ones(n,1),u)-kron(u,ones(n,1));
w=epank(du./h)./h;
W=spdiags(w,0,n^2,n^2);

Yt=kron(ones(n,1),Y);
b=[zeros(q,2*n*p),eye(q)]*((omega'*W*omega)\(omega'*W*Yt));
end


function y = epank( t )
%Epanechnikov kernel
y=0.75*subplus(1-t.^2);
end
