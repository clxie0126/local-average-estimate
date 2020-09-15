function p = T3test(data,I,J)
% GLR2 test
% H0:y=c vs H1:y=f(x)
% data=[u x y]

[n, m]=size(data);
if m-2>I||m-2==I
    error('I should be larger than the number of coefficeint funcion')
end

k=floor(n/I);
n_new=k*I;
data=data(1:n_new,:);
data_sort=sortrows(data);


x=data_sort(:,2:m-1);
y=data_sort(:,m);

X=x(1:I,:);
for i=2:k
    X=blkdiag(X,x(((i-1)*I+1):i*I,:));
end
[l, p]=size(X);
[Q, ~]=qr(X);
v=Q'*y;
u=v(p+1:l);
RSS1=u'*u;

x_j=x(:,J);
x(:,J)=[];
X_dot=x(1:I,:);
for i=2:k
    X_dot=blkdiag(X_dot,x(((i-1)*I+1):i*I,:));
end
D=[X_dot,x_j];

[l, p]=size(D);
[Q, ~]=qr(D);
v=Q'*y;
u=v(p+1:l);
RSS0=u'*u;

T=n_new/2*(RSS0-RSS1)/RSS1;

r=2*(I-m+2)/I;
mu=n_new/I;
p=1-chi2cdf(r*T,mu);

end
