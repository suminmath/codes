%Alg 3.3
function [t1,n1] = alg3_3(x0,u0,D,A,b,p,Eps1,beta,theta,gamma)
t1 = cputime;
n1 = 1;
err1 = 1;
x1n = x0;
u1n = u0;
xn = 2*x1n;
un = 2*u1n;
while err1 > Eps1
   yn = 2*xn-x1n;
   rhon = 2*un-u1n;
   wn = xn+theta*(xn-x1n);
   sigman = un+theta*(un-u1n);
   zn = xn+beta*(xn-x1n);
   pn = wn-gamma*(D'*rhon+A'*(A*zn-b));
   xn1 = max(0,min(pn,1));
   gDyn = gamma*D*yn;
   for i = 1:p
       un1(i) = max(0,sigman(i)+gDyn(i));
   end
   err1 = norm([(xn1-xn)',(un1'-un)'])/norm([xn',un']);
   x1n = xn;
   xn = xn1;
   u1n = un;
   un = un1';
   E1(n1) = err1;
   it1(n1) = n1;
   n1 = n1+1;
end
n1
t1 = cputime-t1
plot(it1,E1,'b','LineWidth',1);
set(gca,'Yscale','log');
xlabel({'\vspace{-1.0mm}';'$n$ '}, 'Interpreter', 'latex');
ylabel({'$E_n$'}, 'Interpreter', 'latex');
legend('Alg 3.3')
hold on


