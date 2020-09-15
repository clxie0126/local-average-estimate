function b = fan(data,h)
%   Fan and Huang(2005) for y=x'a(u)+z'b+e
%   express a(u) in b, then solve b
%   data=[u x z y]
%   local linear fit for coefficient functions
%  

q=1; %the dimension of z
[n, m]=size(data);

u=data(:,1);
x=data(:,2:m-1-q);
z=data(:,m-q:m-1);
Y=data(:,m);

S=zeros(n);
for i=1:n
    U=diag(u-u(i))/h;
    D=[x,U*x];
    W=diag(epank((u-u(i))./h)./h);
    S(i,:)=[x(i,:),zeros(1,m-q-2)]*((D'*W*D)\(D'*W));
end
b=(z'*(eye(n)-S)'*(eye(n)-S)*z)\(z'*(eye(n)-S)'*(eye(n)-S)*Y);
end

function y = epank( t )
%Epanechnikov kernel
y=0.75*subplus(1-t.^2);
end
