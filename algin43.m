%alg (3.9) in [43]
function [t3,n3] = algin43(x0,u0,D,A,b,p,Eps1,gamma)
t3 = cputime;
n3 = 1;
err3 = 1;
x1n = x0;
u1n = u0;
xn = 2*x1n;
un = 2*u1n;
Dx1n = gamma*D*x1n;
Dxn = gamma*D*xn;
while err3 > Eps1
  pn = xn-gamma*(D'*un+A'*(A*xn-b));
  xn1 = max(0,min(pn,1))-gamma*(D'*(un-u1n));
  Dxn1 = gamma*D*xn1;
  for i = 1:p
      un1(i) = max(0,un(i)+Dxn(i))-(Dxn(i)-Dx1n(i));
  end
  err3 = norm([(xn1-xn)',(un1'-un)'])/norm([xn',un']);
  Dx1n = Dxn;
  Dxn = Dxn1；
  u1n = un;
  xn = xn1;
  un = un1';
  E3(n3) = err3;
  it3(n3) = n3;
  n3 = n3+1;
end
n3
t3 = cputime-t3
plot(it3,E3,'Marker','o','MarkerIndices',1:30:length(it3))
set(gca,'Yscale','log');
xlabel({'\vspace{-1.0mm}';'$n$ '}, 'Interpreter', 'latex');
ylabel({'$E_n$'}, 'Interpreter', 'latex');
legend('Alg 3.3', 'alg (2.12) in [7]', 'alg (3.9) in [43]')
hold on
