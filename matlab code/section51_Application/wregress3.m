function sigma2_hat=wregress3(data,a,b)

% refined variance estimator for varyingcoefficient model y=a(u)*x+e
% data=[u,x,y]
% a is the minimum of group size I??
% b is the maximum of group size I??
% so m=b-a+1 is the number of regressions

[n, p]=size(data);
m=b-a+1;

s=zeros(m,1);

d1=zeros(m,1);
d2=zeros(m,1);

for i=1:m
    I=a+i-1;
    k=floor(n/I);
    w=sqrt(k*(I-p+2));
     
    s(i)=average3(data,I)*w;
    d1(i)=w;
    d2(i)=w*(I+1)*(I-1)/((I*k)^2);
end



D=[d1, d2];
beta=regress(s,D);
sigma2_hat=beta(1);
end


function [ sigma ] = average3( data,I )
%local average for varying coefficient model y=a(u)*x+e
% data=[u x y]
[n, m]=size(data);
if m-2>I||m-2==I
    error('I should be larger than the number of coefficeint funcion')
end

k=floor(n/I);
n_new=k*I;
data=data(1:n_new,:);

data_sort=sortrows(data);

y=data_sort(:,m);

X=data_sort(1:I,2:m-1);
for i=2:k
    X=blkdiag(X,data_sort((i-1)*I+1:i*I,2:m-1));
end

[l, p]=size(X);
[Q, ~]=qr(X);
v=Q'*y;
u=v(p+1:l);%QR decomposition, R is upper triangular
sigma=(u'*u)/(n_new-k*(m-2));
end
