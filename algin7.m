%alg (2.12) in [7]
function [t2,n2] = algin7(x0,u0,D,A,b,p,Eps1,gamma)
t2 = cputime;
n2 = 1;
err2 = 1;
un = u0;
xn = x0;
while err2 > Eps1
    pn = xn-gamma*(A'*(A*xn-b)+D'*un);
    yn = max(0,min(pn,1));
    gDxn = gamma*D*xn;
    gDyn = gamma*D*yn;
    m = 0;
    for i = 1:p
        etan(i) = max(0,un(i)+gDxn(i));
        un1(i) = max(0,etan(i)-(gDxn(i)-gDyn(i)));
        m = m+(un(i)-etan(i))*D(i,:);%m¹æÄ£ÊÇ£¨1£¬N£©
    end
    rn = yn+gamma*m';
    xn1 = max(0,min(rn,1));
    err2 = norm([(xn1-xn)',(un1'-un)'])/norm([xn',un']);
    un = un1';
    xn = xn1;
    E2(n2) = err2;
    it2(n2) = n2;
    n2 = n2+1;
end
n2
t2 = cputime-t2
plot(it2,E2,'Marker','*','MarkerIndices',1:30:length(it2))
set(gca,'Yscale','log');
xlabel({'\vspace{-1.0mm}';'$n$ '}, 'Interpreter', 'latex');
ylabel({'$E_n$'}, 'Interpreter', 'latex');
legend('Alg 3.3', 'alg (2.12) in [7]')
hold on
