figure
A = load('Graphics.txt');
A = sort(A);
d1 = A(:, 1);
d2 = A(:, 2);
d2 = log10(d2);
d3 = A(:, 3);
d3 = log10(d3);
plot(d1, d2)
grid on
xlabel('Cond(A)')
%ylabel('log10(||x* - x||)')

hold on
%figure
plot(d1, d3)
grid on
xlabel('Cond(A)')
%ylabel('log10(||Ax - b||)')
legend('||x* - x||', '||Ax - b||')