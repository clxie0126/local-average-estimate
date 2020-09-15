function H = T1test(x,y,h)
% T1 based test
% H0:y=c vs H1:y=f(x)
% The significance level is 0.05

if ~isequal(size(x),size(y))
    error(message('MATLAB:polyfit:XYSizeMismatch'))
end

n=length(x);

c_hat=mean(y);
e=y-c_hat;

S1=0;
S2=0;
for i=1:n
    d=epank((x(i)*ones(n,1)-x)/h).*e;
    d(i)=[];
    S1=S1+sum(d)*e(i);
    S2=S2+sum(d.^2)*e(i)^2;
end
T=S1/sqrt(2*S2);

if T>1.645
   H=1;
else
   H=0;
end
end


function y = epank( t )
%Epanechnikov kernel
y=0.75*subplus(1-t.^2);
end

