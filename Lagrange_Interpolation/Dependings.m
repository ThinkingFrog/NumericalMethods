%Straight on chebyshev and random grid
A = load('output1.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
grid on
hold on
A = load('output3.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
xlabel('Num of points')
ylabel('Max difference')
title('Straight function')
legend('Chebyshev grid', 'Random grid')

%Angled on chebyshev and random grid
figure
A = load('output2.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
grid on
hold on
A = load('output4.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
xlabel('Num of points')
ylabel('Max difference')
title('Angled function at point 0')
legend('Chebyshev grid', 'Random grid')

%Straight and angled on chebyshev grid
figure
A = load('output1.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
grid on
hold on
A = load('output2.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
xlabel('Num of points')
ylabel('Max difference')
title('Cheb grid')
legend('No angle', 'Angle at point 0')

%Straight and angled on random grid
figure
A = load('output3.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
grid on
hold on
A = load('output4.txt');
pnts = A(:, 1);
mstk = A(:, 2);
semilogy(pnts, mstk)
xlabel('Num of points')
ylabel('Max difference')
title('Random grid')
legend('No angle', 'Angle at point 0')