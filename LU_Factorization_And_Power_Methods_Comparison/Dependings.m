figure
A = load('output1.txt');
A = sort(A);
cond = A(:, 1);
precStr = A(:, 2);
precIt = A(:, 3);
precItLess = A(:, 4);
precItMore = A(:, 5);
grid on
hold on
xlabel('cond(A)')
ylabel('||x* - x||')
plot(cond, precStr)
plot(cond, precIt)
plot(cond, precItLess)
plot(cond, precItMore)
legend('||x* - x||(LU)', '||x* - x||(iter)', '||x* - x||(iter) / 2', '||x* - x||(iter) * 2')

figure
A = load('output2.txt');
eps = A(:, 1);
str1 = A(:, 2);
str2 = A(:, 3);
it1 = A(:, 4);
it2 = A(:, 5);
semilogx(eps, str2)
hold on
semilogx(eps, str2 / 2)
semilogx(eps, str2 * 2)
semilogx(eps, it2)
grid on
legend('N(LU)', 'N(LU) / 2', 'N(LU) * 2', 'N(iter)')
xlabel('epsilon')
ylabel('N')