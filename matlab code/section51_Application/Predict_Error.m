%% the dataset is included in UsedData.mat
%% the first column is the vector of time t, which is the variable U in this case
%% the last column is the vector of response variable.
%% the column 2 is an one vector, which corresponds to the intercept. 
%% the column 3,4,5 correpond to X_2, X_3, X_4.

load('UsedData.mat');
SortData=sortrows(UsedData); u = SortData(:,1); y = SortData(:,6);

%% the paramter I is 10
%% the bandwidth h is taken to be 30% of the interval length of U.
n = length(u); n_train = 300; I = 10; h = 0.3*300; 
k=floor(n_train/I); n_new=k*I;

%% the sample size of the train data is taken to be 500. 
%% the local average estimate procedure
PredictValue = zeros((n - n_train),1);
for s=1:(n - n_train)
    
x_train = SortData( s:(s+n_train-1) ,2:5); 
u_train = SortData(s:(s+n_train-1),1);
y_train = SortData(s:(s+n_train-1),6);

u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u_train((I*i-I+1):I*i));
end
X=x_train(1:I,:);
for i=2:k
    X=blkdiag(X,x_train((I*i-I+1):i*I,:));
end
b=regress(y_train,X);  %estimated a(u_star)
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
clear P; clear W; clear x_train; clear y_train; clear u_train;

%% predict value
PredictValue(s) = SortData(s+n_train,2:5)*au(k,:)';

end
TMAE = sum(abs( y((n_train+1):730) - PredictValue))/(730 - n_train);
NMAE = sum(abs( y((n_train+1):730) - y(n_train:729)))/(730 - n_train);
TMSE = norm( y((n_train+1):730) - PredictValue)^2/(730 - n_train);
NMSE = norm( y((n_train+1):730) - y(n_train:729))^2/(730 - n_train);

Tresult = [PredictValue, y((n_train+1):730) - PredictValue ];


%% predict error of reduced model

PredictValue = zeros(n - n_train,1);
for s=1:(n - n_train)
    
x_train = SortData( s:(s+n_train-1) ,[2,5]); 
u_train = SortData(s:(s+n_train-1),1);
y_train = SortData(s:(s+n_train-1),6);

u_star=zeros(k,1);
for i=1:k
    u_star(i)=mean(u_train((I*i-I+1):I*i));
end
X=x_train(1:I,:);
for i=2:k
    X=blkdiag(X,x_train((I*i-I+1):i*I,:));
end
b=regress(y_train,X);  %estimated a(u_star)
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

%% predict value
PredictValue(s) = SortData(s+n_train,[2,5])*au(k,:)';

end

RMAE = sum(abs( y((n_train+1):730) - PredictValue))/(730 - n_train);
RMSE = norm( y((n_train+1):730) - PredictValue)^2/(730 - n_train);
Rresult = [PredictValue, y((n_train+1):730) - PredictValue ];

%% the mean absolute error of the total(full) model, reduced model and null model

TMAE
RMAE
NMAE

