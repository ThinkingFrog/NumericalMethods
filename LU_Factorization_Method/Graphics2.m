A = load('Graphics2.txt');
d1 = A(:, 1);
%d1 = log10(d1);
d2 = A(:, 2);
plot(d2, d1, '.')
grid on
xlabel('vector error in %')
ylabel('||x* - x|| / ||x*||')

hold on
x = 0 : 0.01 : 5;
y = x * 1030;
%y = log10(y);
plot(x, y)
legend('||x* - x|| / ||x*||', '||b* - b|| / ||b*|| * cond(A)')