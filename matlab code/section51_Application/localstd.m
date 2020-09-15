function [u_star, u_std] = localstd(data,I,d)
%   upperbound of local average for y=x'a(u)+e
%   average u
%   data=[u x y]
%   now p means the polynomial order in second step ployfit
%  NOTICE!!! the bandwith constant h or nearest d


j=2;

%estimate residual variance by local average
sigma2_hat=wregress3(data,5,14);


[n, m]=size(data);
if m-2>I||m-2==I
    error('I should be larger than the number of coefficeint funcion')
end

k=floor(n/I);
n_new=k*I;
data=data(1:n_new,:);

data_sort=sortrows(data);

u=data_sort(:,1);
u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u(I*(i-1)+1:I*i));
end

sigma2_u=zeros(k,1);
for i=1:k
    X=data_sort((i-1)*I+1:i*I,2:m-1);
    fai=I*inv(X'*X);
    sigma2_u(i)=unit(j,m-2)*fai*unit(j,m-2)';
end


h=d*729;
u_std1=sqrt(1.25*sigma2_hat*sigma2_u./h);

%sigma2_m=zeros(k,1);

%for i=1:k
    %U=[ones(k,1), u_star-u_star(i),(u_star-u_star(i)).^2,(u_star-u_star(i)).^3 ];
    %W=diag(epank((u_star-u_star(i))./h)./h);
    %sigma2_m(i)=unit(1,4)*(inv(U'*W*U)*(U'*W*W*U)/(U'*W*U))*unit(1,4)';
%end

%u_std1=sqrt(sigma2_m.*sigma2_u.*sigma2_hat);

w = weight( u_star,h );
P = local_poly( u_star,u_std1,3,w );
u_std=P(:,4);

end

