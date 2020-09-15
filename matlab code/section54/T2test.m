function H = T2test(x,y)
% GLR test
% H0:y=c vs H1:y=f(x)
% f(x) using local linear estimator
% The significant level is 0.05

if ~isequal(size(x),size(y))
    error(message('MATLAB:polyfit:XYSizeMismatch'))
end

n=length(x);h=n^(-1/5);


c_hat=mean(y);
RSS0=sum((y-c_hat).^2);

W = zeros(n,n);
for k = 1:n
    Xk = x(k)*ones(size(x));
    W(k,:) = (epank((x-Xk)./h)./h)';
end
P = zeros(n,2);
for k = 1:n
    Xk = x(k)*ones(size(x));
    x_k = x - Xk;
    P(k,:)= wpolyfit(x_k,y,1,W(k,:));
end
fx=P(:,2);
RSS1=sum((y-fx).^2);

T=(n/2)*log(RSS0/RSS1);

%r=2.1153;
%omega=1;
mu=2.1153*0.45*1/h;
chi2=chi2inv(0.95,mu);

if 2.1153*T>chi2
   H=1;
else
   H=0;
end

end

function y = epank( t )
%Epanechnikov kernel
y=0.75*subplus(1-t.^2);
end