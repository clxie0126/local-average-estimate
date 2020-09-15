%% the dataset is included in UsedData.mat
%% the first column is the vector of time t, which is the variable U in this case
%% the last column is the vector of response variable.
%% the column 2 is an one vector, which corresponds to the intercept. 
%% the column 3,4,5 correpond to X_2, X_3, X_4.
load('UsedData.mat');
SortData=sortrows(UsedData);
u = SortData(:,1);
y = SortData(:,6);

%% the paramter I is 10
%% the bandwidth h is taken to be 30% of the interval length of U.
n=length(u); I =10; h = 0.3*729;
k=floor(n/I); n_new=k*I;

%% the local average estimate procedure
u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u((I*i-I+1):I*i));
end
X=SortData(1:I,2:5);
for i=2:k
    X=blkdiag(X,SortData((I*i-I+1):i*I,2:5));
end
b=regress(y,X);  %estimated a(u_star)
B=reshape(b,4,k);
clear b;

%% smooth step after local average estimation
W = 0.75*subplus( 1-(( u_star*ones(1,k) - ones(k,1)*u_star' )/h).^2);
au = zeros(k,4);
for i=1:4
    aui_star=B(i,:);
    P = local_poly( u_star,aui_star, 3,W );
    au(:,i)=P(:,3+1);
end
clear P;

%% estimate residual variance by local average: 'sigma2_hat'
%% estimate \Gamma
%% estimate the variance of \hat_au
sigma2_hat=wregress3(SortData,5,14);
sigma2_u=zeros(k,4);
for j=1:4
unit_j4 = zeros(1,4);
unit_j4(j) = 1;
for i=1:k
    xi=SortData((i-1)*I+1:i*I,2:5);
    fai=I*inv(xi'*xi);
    sigma2_u(i,j)=unit_j4*fai*unit_j4';
end
end

%sigma2_u = diag(I*eys(4*k)/(X'*X));
%sigma2_u = reshape(sigma2_u, 4,k);
%sigma2_u = sigma2_u';
u_std0=(1.25*sigma2_hat*sigma2_u./h).^(1/2);

u_std = zeros(k,4);
for j=1:4
P = local_poly( u_star,u_std0(:,j),3,W );
u_std(:,j)=P(:,4);
end

%% the pointwise confidence interval of \hat_au:(\hat_au-1.96std,\hat_au+1.96std)
%% the smoothed pointwise confidence interval: smooth the upper bound and lower bound respectively
UpperCI_au = au + 1.96*u_std;
LowerCI_au = au - 1.96*u_std;

SmUpperCI_au = zeros(k,4);
for j=1:4
P = local_poly( u_star,UpperCI_au(:,j),3,W);
SmUpperCI_au(:,j)=P(:,4);
end

SmLowerCI_au = zeros(k,4);
for j=1:4
P = local_poly( u_star,LowerCI_au(:,j),3,W);
SmLowerCI_au(:,j)=P(:,4);
end

clear P;

%% plot the estimate of au and its smoothed pointwise confidence bounded 
ax1 = subplot(2,2,1);
plot(ax1, au(:,1))
hold on
plot(ax1, UpperCI_au(:,1))
hold on
plot(ax1, LowerCI_au(:,1))
xlabel(ax1,'t')
ax2 = subplot(2,2,2);
plot(ax2, au(:,2))
hold on
plot(ax2, UpperCI_au(:,2))
hold on
plot(ax2, LowerCI_au(:,2))
xlabel(ax2,'t')
ax3 = subplot(2,2,3);
plot(ax3, au(:,3))
hold on
plot(ax3, UpperCI_au(:,3))
hold on
plot(ax3, LowerCI_au(:,3))
xlabel(ax3,'t')
ax4 = subplot(2,2,4);
plot(ax4, au(:,4))
hold on
plot(ax4, UpperCI_au(:,4))
hold on
plot(ax4, LowerCI_au(:,4))
xlabel(ax4,'t')

%% test procedure used T_3
%% the null hypothesis is au =0

[l, p]=size(X); [Q, ~]=qr(X);
v=Q'*y; u=v(p+1:l);
RSS1=u'*u;
clear p; clear Q; clear l; clear v; clear u;

r=2*(I-4)/I;
mu=n_new/I;
Pvalue = zeros(1,4);
for j=1:4
    D=SortData(:,2:5);
    D(:,j)=[];
    X_dot=D(1:I,:);
    for i=2:k
        X_dot=blkdiag(X_dot,D((I*i-I+1):i*I,:));
    end
    [l, p]=size(X_dot); [Q, ~]=qr(X_dot);
    v=Q'*y; u=v(p+1:l); RSS0=u'*u;
    T=n_new/2*(RSS0-RSS1)/RSS1;
    Pvalue(j)=1-chi2cdf(r*T,mu);
    clear D; clear p; clear Q; clear l; clear v; clear u;
end

%% the test procedure used T_3
%% the null hypothesis is 'au = c'

r=2*(I-4)/I;
mu=n_new/I;
Pvalue = zeros(1,4);
for j=1:4
    D=SortData(:,2:5);
    d_j = D(:,j);
    D(:,j)=[];
    X_dot=D(1:I,:);
    for i=2:k
        X_dot=blkdiag(X_dot,D((I*i-I+1):i*I,:));
    end
    Z_dot = [X_dot, d_j];
    [l, p]=size(Z_dot); [Q, ~]=qr(Z_dot);
    v=Q'*y; u=v(p+1:l); RSS0=u'*u;
    T=n_new/2*(RSS0-RSS1)/RSS1;
    Pvalue(j)=1-chi2cdf(r*T,mu);
    clear D; clear p; clear Q; clear l; clear v; clear u;
end




