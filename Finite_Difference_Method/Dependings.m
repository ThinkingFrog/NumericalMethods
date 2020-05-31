A = load('output1.txt');
h = A(:, 1);
mstk = A(:, 2);
loglog(h, mstk)
grid on
xlabel('h')
ylabel('Max mistake')
title('|y*-y|(h)')

figure
A = load('output2.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)
grid on
xlabel('delta')
ylabel('Max mistake')
title('Left part delta')

figure
A = load('output3.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)
grid on
xlabel('delta')
ylabel('Max mistake')
title('Right part delta')