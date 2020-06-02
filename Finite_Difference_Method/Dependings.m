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

hold on
A = load('output3.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)

legend('Left part delta', 'Right part delta')
title('|y*-y|(delta)')