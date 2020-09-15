function b = zhang(data,h)
%   Zhang(2002) for y=x'a(u)+z'b+e
%   treat b as function, then take average
%   data=[u x z y]
%   local linear fit for coefficient functions
%  

q=1; %the dimension of z
J=1; %the Jth constant coefficient 

[n, m]=size(data);

u=data(:,1);
x=data(:,2:m-1-q);
z=data(:,m-q:m-1);
Y=data(:,m);

bu=zeros(n,1);
l=(m-q-2)*2+q;
for i=1:n
U=diag(u-u(i));
V=[z,x,U*x];
W=diag(epank((u-u(i))./h)./h);
bu(i)=[zeros(1,J-1),1,zeros(1,l-J)]*((V'*W*V)\(V'*W*Y));
end
b=mean(bu);
end

function y = epank( t )
%Epanechnikov kernel
y=0.75*subplus(1-t.^2);
end

