%% Example in Section 5.3
%% Fan and Huang(2005) for y=x'a(u)+z'b+e

n=500;
sigma=[1,0,0;0,1,0;0,0,1];

Estb = zeros(100,1);
for rep = 1:100
    U=rand(n,1);
    X=mvnrnd([0,0,0],sigma,n);
    Y_mean=diag(X*[sin(2*pi*U),cos(2*pi*U),ones(n,1)]'); %Example 4
    %Y_mean=diag(X*[sin(2*pi*U),3.5*(exp(-(4*U-1).^2)+exp(-(4*U-3).^2))-1.5,ones(n,1)]'); %Example 5
    %Y_mean=diag(X*[sin(6*pi*U),sin(2*pi*U),ones(n,1)]'); %Example 6
    sd2=0.2*var(Y_mean);
    e=normrnd(0,sqrt(sd2),n,1);
    Y=Y_mean+e;
    data=[U,X,Y];
    Estb(rep) = fan(data,0.5);
end

mean(Estb)
std(Estb)