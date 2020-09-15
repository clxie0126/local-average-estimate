function [u_star, au] = locala_varying(data, Sorder, I, h)
%   local average for varying coefficient model: y=x'a(u)+e
%   data=[u x y], where u is a n-length vector, x is a n by p matrix and y
%   is a n-length vector
%  The parameter Sorder(smooth order) is the polynomial order in second step ployfit


[n, m]=size(data); p = m-2;
if p>I||p==I
    error('I should be larger than the number of coefficeint funcion')
end

k=floor(n/I); nn=k*I;
data=data(1:nn,:);
data_sort=sortrows(data);

u=data_sort(:,1);
y=data_sort(:,m);

u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u((I*i-I+1):I*i));
end
X=data_sort(1:I,2:m-1);
for i=2:k
    X=blkdiag(X,data_sort((i*I-I+1):i*I,2:m-1));
end

b=regress(y,X); au_star=reshape(b,p,k); au_star = au_star';
W = 0.75*subplus( 1-(( u_star*ones(1,k) - ones(k,1)*u_star' )/h).^2);
au = zeros(k, p);
for j=1:p
    P = local_poly( u_star, au_star(:,j), Sorder, W );
    au(:,j)=P(:,Sorder+1);
end


end

