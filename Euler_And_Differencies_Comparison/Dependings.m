%1st plot
A = load('output1_1.txt');
h = A(:, 1);
mstk = A(:, 2);
loglog(h, mstk)

hold on
grid on

A = load('output1_2.txt');
h = A(:, 1);
mstk = A(:, 2);
loglog(h, mstk)

h = [1e-4, 1e-0];
mstk = h .^ 2;
loglog(h, mstk)

xlabel('h')
ylabel('|y*-y|')
legend('Euler-Cauchy', 'MKP', 'h^2')
title('|y*-y|(h)')

%2nd plot
figure
A = load('output2_1.txt');
h = A(:, 1);
coord = A(:, 2);
semilogx(h, coord)

hold on
grid on

A = load('output2_2.txt');
h = A(:, 1);
coord = A(:, 2);
semilogx(h, coord)

xlabel('h')
ylabel('x')
legend('Euler-Cauchy', 'MKP')
title('x(h)')

%3rd plot
figure
A = load('output3_1.txt');
h = A(:, 1);
V = A(:, 2);
loglog(h, V)

hold on
grid on

A = load('output3_2.txt');
h = A(:, 1);
V = A(:, 2);
loglog(h, V)

xlabel('h')
ylabel('V')
legend('Euler-Cauchy', 'MKP')
title('V(h)')

%4th plot
figure
A = load('output4_1.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)

hold on
grid on

A = load('output4_2.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)

xlabel('delta1')
ylabel('|y*-y|')
legend('Euler-Cauchy', 'MKP')
title('|y*-y|(delta1)')

%5th plot
figure
A = load('output5_1.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)

hold on
grid on

A = load('output5_2.txt');
delta = A(:, 1);
mstk = A(:, 2);
loglog(delta, mstk)

xlabel('delta2')
ylabel('|y*-y|')
legend('Euler-Cauchy', 'MKP')
title('|y*-y|(delta2)')