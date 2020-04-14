A = load('output1.txt');
mistake = A(:, 1);
eps = A(:, 2);
max = A(:, 3);
min = A(:, 4);
loglog(eps, mistake)
grid on
hold on
loglog(eps, min)
hold on
loglog(eps, max, 'g')
legend('my', 'min', 'max')
xlabel('eps')
%ylabel('|I* - I|')
title('|I* - I|(eps)')

figure
A = load('output2.txt');
L = A(:, 1);
eps = A(:, 2);
semilogx(eps, L)
grid on
xlabel('eps')
ylabel('L')
title('M = 2^L (eps)')