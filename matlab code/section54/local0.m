function [u_star, au_star] = local0(data,I,j)
%   local average for y=x'a(u)+e
%   average u
%   data=[u x y]
%   primary point estiamtors for jth coefficient


[n, m]=size(data);
if m-2>I||m-2==I
    error('I should be larger than the number of coefficeint funcion')
end

k=floor(n/I); nn=k*I;
data=data(1:nn,:);
data_sort=sortrows(data);

u=data_sort(:,1);
y=data_sort(:,m);

u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u(I*(i-1)+1:I*i));
end
X=data_sort(1:I,2:m-1);
for i=2:k
    X=blkdiag(X,data_sort((i-1)*I+1:i*I,2:m-1));
end

b=regress(y,X);  %estimated a(u_star)
B=reshape(b,m-2,k);
au_star=B(j,:);
au_star=au_star';

end

