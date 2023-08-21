clc
clear
Eps1 = 1e-6;
sum1 = 0;
sum2 = 0;
sum3 = 0;
sum4 = 0;
sum11 = 0;
sum22 = 0;
sum33 = 0;
sum44 = 0;
m = 2500;
N = 2*m;
p = 600;
delta = 4;
theta = 0.01;
beta = 0.01;
epsilon3 = 2.7;%(2,3) 
epsilon1 = 0.99*(0.5-1/epsilon3);
Ia = 10; 
for f = 1: Ia
    A = 10*rand(m,N);
    b = 5*rand(m,1);
    D = 8*rand(p,N);
    x0 = rand(N,1);
    u0 = rand(p,1);
    l1 = norm(D);
    l2 = norm(A)^2;
    xi = (delta/l2)/(1+(1+16*(1/l2)^2*l1^2)^0.5);
    eta = 0.01*xi;%(0,0.5*xi)
    nu = l1^2*eta*(2*xi-eta)/(1+(xi-eta)*l1)^2;
    alphan = 0.99*min(((1+8*nu)^0.5-1-2*nu)/(2*(1-nu)),1);
    box1 = 0.99*(1-4*theta)/(3*l2*(2+2*beta+beta^2)+(1+2^0.5)*l1);
    box2 = 0.99*xi;%[eta,xi-eta]
    epsilon2 = 0.6/l2;
    up1 = 0.99*(2/l2-epsilon2)*epsilon1;
    up2 = (3-epsilon3)*epsilon2;
    up3 = 0.99*(0.5-epsilon1-1/epsilon3)/l1;
    box3 = min([up1,up2,up3]);
    gamma = min([box1,box2,box3]);
    [t1,n1] = alg3_3(x0,u0,D,A,b,p,Eps1,beta,theta,gamma);
    [t2,n2] = algin7(x0,u0,D,A,b,p,Eps1,gamma);
    [t3,n3] = algin43(x0,u0,D,A,b,p,Eps1,gamma);
    [t4,n4] = algin46(x0,u0,D,A,b,p,Eps1,alphan,gamma);
    sum1 = sum1+t1;
    sum11 = sum11+n1;
    sum2 = sum2+t2;
    sum22 = sum22+n2;
    sum3 = sum3+t3;
    sum33 = sum33+n3;
    sum4 = sum4+t4;
    sum44 = sum44+n4;
 end
ave1 = sum1/Ia
ave11 = sum11/Ia
ave2 = sum2/Ia
ave22 = sum22/Ia
ave3 = sum3/Ia
ave33 = sum33/Ia
ave4 = sum4/Ia
ave44 = sum44/Ia

