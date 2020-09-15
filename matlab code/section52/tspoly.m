function au = tspoly( data,j,h2 )
%two step local poly,  for varying coefficeint model
% y = X'a(u)+e
%Fan and Zhang 1999

%data=[u x y], j is for a_j(u)
%au = [a_j(u1), a_j(u2), ..., a_j(un)]; R^n vector

[n,m]=size(data);

if mod(j,1)
   error('j should be an integer')
elseif j<1||j>m-2
   error('j is out of range')
end

X=data(:,2:m-1);
Y=data(:,m);
u=data(:,1);

%initial estimator,local linear with h0=0.05
au0 = initial(data,0.05);

%second step, local cubic for a_p(u) with h2
au0(j,:)=[];
X0=X;
X0(:,j)=[];
V=diag(X0*au0);

x=X(:,j);
au=zeros(n,1);
for i=1:n
U=diag(u-u(i));
X2=[x,U*x,U^2*x,U^3*x];
W=diag(epank((u-u(i))./h2)./h2);
w=sqrt(W);
y_delta=w*(Y-V);
x_delta=w*X2;
[Q,R] = qr(x_delta,0);
p = R\(Q'*y_delta);
au(i)=[1,0,0,0]*p;
end

end

function au = initial(data,h)
%one step local linear,  for varying coefficeint model
% y = X'a(u)+e
%Fan and Zhang 2008
%data=[u x y]
% au = [a_1(u)
%       a_2(u)
%       
%       a_p(u)];   p*n matrix



if nargin < 2
    h='default';
end

[n,m]=size(data);

X=data(:,2:m-1);
Y=data(:,m);
u=data(:,1);

switch h
    case 'default'
        h=band_width(u,0.3);
end

au=zeros(m-2,n);

for i=1:n
U=diag(u-u(i));
gama=[X,U*X];
W=diag(epank((u-u(i))./h)./h);
au(:,i)=[eye(m-2),zeros(m-2)]*((gama'*W*gama)\(gama'*W*Y));
end

end

function y = epank( t )
y=0.75*subplus(1-t.^2);
end

