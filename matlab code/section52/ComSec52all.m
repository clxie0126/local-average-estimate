%% Section 5.2
%% Example 1
n = 250; I = 10; k = n/I;
MISE_h = zeros(5,3);
for j=1:5
    h= j*0.2;
    MISE = zeros(10,3);
for rep =1:10
    U=rand(n,1);
    X=mvnrnd([0,0],[1,sqrt(1/2);sqrt(1/2),1],n);
    Y_mean=diag(X*[sin(60*U),4*U.*(1-U)]');     %EXAMPLE 1
    sigma2=0.2*var(Y_mean);
    e=normrnd(0,sqrt(sigma2),n,1);
    Y=Y_mean+e;
    data=[U,X,Y];
    [u_star, au_lav] = locala_varying( data, 3, I, h);
    AU= 4*u_star.*(1-u_star); 
    MISE(rep,1)=norm(au_lav(:,2)-AU)^2/k;
    au_one = ospoly(data,2,h);
    au_two = tspoly(data,2,h);
    AU=4*U.*(1-U);  
    MISE(rep,2)=norm(au_one-AU)^2/n;
    MISE(rep,3)=norm(au_two-AU)^2/n;
end
MISE_h(j,:)=mean(MISE);
j
end

%% Example 2

n = 500; I = 10; k = n/I;
MISE_h = zeros(10,3);
for j=1:10
    h= j/10;
    MISE = zeros(10,3);
for rep =1:10
    U=rand(n,1);
    X=mvnrnd([0,0],[1,sqrt(1/2);sqrt(1/2),1],n);
    Y_mean=diag(X*[sin(6*pi*U),sin(2*pi*U)]');    %EXAMPLE 2
    sigma2=0.2*var(Y_mean);
    e=normrnd(0,sqrt(sigma2),n,1);
    Y=Y_mean+e;
    data=[U,X,Y];
    [u_star, au_lav] = locala_varying( data, 3, I, h);
    AU=sin(2*pi*u_star);
    MISE(rep,1)=norm(au_lav(:,2)-AU)^2/k;
    au_one = ospoly(data,2,h);
    au_two = tspoly(data,2,h);
    AU=sin(2*pi*U);  
    MISE(rep,2)=norm(au_one-AU)^2/n;
    MISE(rep,3)=norm(au_two-AU)^2/n;
end
MISE_h(j,:)=mean(MISE);
j
end

%% Example 3

n = 500; I = 10; k = n/I;
MISE_h = zeros(10,3);
for j=1:10
    h= j/10;
    MISE = zeros(10,3);
for rep =1:10
    U=rand(n,1);
    X=mvnrnd([0,0],[1,sqrt(1/2);sqrt(1/2),1],n);
    Y_mean=diag(X*[sin(8*pi*(U-0.5)),(3.5*(exp(-(4*U-1).^2)+exp(-(4*U-3).^2))-1.5)]');     %EXAMPLE 3
    sigma2=0.2*var(Y_mean);
    e=normrnd(0,sqrt(sigma2),n,1);
    Y=Y_mean+e;
    data=[U,X,Y];
    [u_star, au_lav] = locala_varying( data, 3, I, h);
    AU=3.5*(exp(-(4*u_star-1).^2)+exp(-(4*u_star-3).^2))-1.5; 
    MISE(rep,1)=norm(au_lav(:,2)-AU)^2/k;
    au_one = ospoly(data,2,h);
    au_two = tspoly(data,2,h);
    AU=3.5*(exp(-(4*U-1).^2)+exp(-(4*U-3).^2))-1.5;   
    MISE(rep,2)=norm(au_one-AU)^2/n;
    MISE(rep,3)=norm(au_two-AU)^2/n;
end
MISE_h(j,:)=mean(MISE);
j
end
