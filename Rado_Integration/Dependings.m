A = load('output1.txt');
eps = A(:, 1);
mstk = A(:, 2);
loglog(eps, mstk)
grid on
hold on
xlabel('eps')
ylabel('mstk')
title('|I*-I|(eps)')

figure
A = load('output2.txt');
eps = A(:, 1);
L = A(:, 2);
semilogx(eps, L)
grid on
hold on
xlabel('eps')
ylabel('L')
title('M = 2^L (eps)')

figure 
A = load('output3.txt');
eps = A(:, 1);
mstk = A(:, 2);
h = A(:, 3);
loglog(eps, mstk)
grid on
hold on
loglog(eps, h)
xlabel('eps')
ylabel('mstk')
legend('|I*-I|', 'h^p')