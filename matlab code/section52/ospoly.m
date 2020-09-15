function au = ospoly(data,j,h)
%one step local polynomial,  for varying coefficeint model
% y = X'a(u)+e
%Fan and Zhang 2008
%data=[u x y]
%j means jth function
%au = [a_j(u1), a_j(u2), ..., a_j(un)]; R^n vector



[n,m]=size(data);

X=data(:,2:m-1);
Y=data(:,m);
u=data(:,1);

x=X(:,j);
X(:,j)=[];
t=2*m-6;

au=zeros(n,1);
for i=1:n
U=diag(u-u(i));
X2=[x,U*x,U^2*x,U^3*x];
X3=[X,U*X];
X1=[X3,X2];
W=diag(epank( ( u-u(i) )./h )./h);
au(i)=[zeros(1,t),1,0,0,0]*((X1'*W*X1)\(X1'*W*Y));
end

end

function y = epank( t )
y=0.75*subplus(1-t.^2);
end

