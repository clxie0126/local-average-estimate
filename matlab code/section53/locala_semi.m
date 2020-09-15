function [b] = locala_semi(data, I)
%   local average for semivarying coefficient model b in y=x'a(u)+z'b+e,
%   estimate b
%   data=[u x z y], where u is a n-length vector, x is a n by p matrix and y
%   is a n-length vector
%   The parameter Sorder(smooth order) is the polynomial order in second step ployfit

[n, m]=size(data);
if m-2>I||m-2==I
    error('I should be larger than the number of coefficeint funcion')
end

k=floor(n/I); nn=k*I;
data=data(1:nn,:);
data_sort=sortrows(data);

u=data_sort(:,1);
z=data_sort(:,(m-1));
y=data_sort(:,m);

u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u((I*i-I+1):I*i));
end
X=data_sort(1:I,2:m-2);
for i=2:k
    X=blkdiag(X,data_sort((i*I-I+1):i*I,2:m-2));
end
D=[z,X];
hat_th=regress(y,D); 
b = hat_th(1);

end

