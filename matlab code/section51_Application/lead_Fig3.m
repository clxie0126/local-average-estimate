%% the dataset is included in UsedData.mat
%% the first column is the vector of time t, which is the variable U in this case
%% the last column is the vector of response variable.
%% the column 2 is an one vector, which corresponds to the intercept. 
%% the column 3,4,5 correpond to X_2, X_3, X_4.
%% in this case, we only consider X_1 and X_4
load('UsedData.mat');
UsedData(:,[3,4]) = [];
SortData=sortrows(UsedData);
u = SortData(:,1);
y = SortData(:,4);

%% the paramter I is 10
%% the bandwidth h is taken to be 30% of the interval length of U.
n=length(u); I =10; h = 0.3*729;
k=floor(n/I); n_new=k*I;

%% the local average estimate procedure
u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u((I*i-I+1):I*i));
end
X=SortData(1:I,2:3);
for i=2:k
    X=blkdiag(X,SortData((I*i-I+1):i*I,2:3));
end
b=regress(y,X);  %estimated a(u_star)
B=reshape(b,2,k);
clear b;

%% smooth step after local average estimation
W = 0.75*subplus( 1-(( u_star*ones(1,k) - ones(k,1)*u_star' )/h).^2);
au = zeros(k,2);
for i=1:2
    aui_star=B(i,:);
    P = local_poly( u_star,aui_star, 3,W );
    au(:,i)=P(:,3+1);
end
clear P;

%% estimate residual variance by local average: 'sigma2_hat'
%% estimate \Gamma
%% estimate the variance of \hat_au
sigma2_hat=wregress3(SortData,5,14);
sigma2_u=zeros(k,2);
for j=1:2
unit_j2 = zeros(1,2);
unit_j2(j) = 1;
for i=1:k
    xi=SortData((i-1)*I+1:i*I,2:3);
    fai=I*inv(xi'*xi);
    sigma2_u(i,j)=unit_j2*fai*unit_j2';
end
end
u_std0=(1.25*sigma2_hat*sigma2_u./h).^(1/2);

u_std = zeros(k,2);
for j=1:2
P = local_poly( u_star,u_std0(:,j),3,W );
u_std(:,j)=P(:,4);
end

%% the pointwise confidence interval of \hat_au:(\hat_au-1.96std,\hat_au+1.96std)
%% the smoothed pointwise confidence interval: smooth the upper bound and lower bound respectively
UpperCI_au = au + 1.96*u_std;
LowerCI_au = au - 1.96*u_std;

SmUpperCI_au = zeros(k,2);
for j=1:2
P = local_poly( u_star,UpperCI_au(:,j),3,W);
SmUpperCI_au(:,j)=P(:,4);
end

SmLowerCI_au = zeros(k,2);
for j=1:2
P = local_poly( u_star,LowerCI_au(:,j),3,W);
SmLowerCI_au(:,j)=P(:,4);
end

clear P;

%% plot the estimate of au and its smoothed pointwise confidence bounded 
%% Problem: Also the confidence bounds is different to those in the paper !!!
%% what is the used bandwidth for the smoothing step?
ax1 = subplot(1,2,1);
plot(ax1, au(:,1))
hold on
plot(ax1, UpperCI_au(:,1))
hold on
plot(ax1, LowerCI_au(:,1))
xlabel(ax1,'t')
ax2 = subplot(1,2,2);
plot(ax2, au(:,2))
hold on
plot(ax2, UpperCI_au(:,2))
hold on
plot(ax2, LowerCI_au(:,2))
xlabel(ax2,'t')
