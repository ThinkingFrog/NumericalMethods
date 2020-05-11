A = load('output1.txt');
eps = A(:, 1);
mstk = A(:, 2);
loglog(eps, mstk)
grid on
xlabel('epsilon')
ylabel('Max mistake')
title('|y*-y|(eps)')

figure
A = load('output2.txt');
eps = A(:, 1);
divs = A(:, 2);
divs = log2(divs)
semilogx(eps, divs)
grid on
xlabel('epsilon')
ylabel('number of iterations')
title('iterations(eps)')

figure
A = load('output3.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)
grid on
xlabel('delta')
ylabel('max mistake')
title('|y*-y|(delta)')