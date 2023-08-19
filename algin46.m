% Alg 1 in [46]
function [t4,n4] = algin46(x0,u0,D,A,b,p,Eps1,alphan,gamma)
t4 = cputime;
n4 = 1;
err4 = 1;
z1n = 1*x0;
u1n = 1*u0;
zn = 2*x0;
un = 2*u0;
 while err4 > Eps1
   wn = zn+alphan*(zn-z1n);
   mn = un+alphan*(un-u1n);%(分块)
   xnn = wn-gamma*(D'*mn+A'*(A*wn-b));
   xn = max(0,min(xnn,1));
   gDxn = gamma*D*xn;
   gDwn = gamma*D*wn;
   m = 0;
   for i = 1:p
       gn(i) = max(0,mn(i)+gDwn(i));
       un1(i) = gn(i)-(gDwn(i)-gDxn(i));
       m = m+(mn(i)-gn(i))*D(i,:);%m规模是（1，N）
   end
   zn1 = xn+gamma*m';
   err4 = norm([(zn1-zn)',(un1'-un)'])/norm([zn',un']);
   z1n = zn;
   zn = zn1;
   u1n = un;
   un = un1';
   E4(n4) = err4;
   it4(n4) = n4;
   n4 = n4+1;
end
n4
t4 = cputime-t4
plot(it4,E4,'g--','LineWidth',1);
set(gca,'Yscale','log');
xlabel({'\vspace{-1.0mm}';'$n$ '}, 'Interpreter', 'latex');
ylabel({'$E_n$'}, 'Interpreter', 'latex');
legend('Alg 3.3', 'alg (2.12) in [7]','alg (3.9) in [43]', 'Alg 1 in [46]')
hold on