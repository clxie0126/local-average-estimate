%% compute emprtical size and power of T2

N=100; n = 400; % n=800; n=1600;
I=10; %I=4; I=5;

RejectTime = 0;

for i=1:N
    U=rand(n,1);
    X=mvnrnd([0,0],[1,sqrt(1/2);sqrt(1/2),1],n);
    Y_mean=diag(X*[sin(60*U),ones(size(U))]');  %null model
    e=normrnd(0,0.557,n,1);
    Y=Y_mean+e;
    data=[U,X,Y];
    p = T3test(data,I,2);
    RejectTime = RejectTime + (p<0.05);
    i
end
a=RejectTime/N