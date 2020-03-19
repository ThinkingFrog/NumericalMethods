figure
A = load('output1_1.txt');
A = sort(A);
lambda = A(:, 1);
eps = A(:, 2);
loglog(eps, lambda)
grid on
xlabel('eps')
ylabel('||lambda* - lambda||')

figure
A = load('output1_2.txt');
A = sort(A);
nevyazka = A(:, 1);
eps = A(:, 2);
loglog(eps, nevyazka)
grid on
xlabel('eps')
ylabel('||Ax - lambda * x||')

figure
A = load('output1_3.txt');
A = sort(A);
x = A(:, 1);
eps = A(:, 2);
semilogx(eps, x)
grid on
xlabel('eps')
ylabel('||x* - x||')

figure
A = load('output2_1.txt');
A = sort(A);
lambda = A(:, 1);
move = A(:, 2);
plot(move, lambda)
grid on
xlabel('delta')
ylabel('||lambda* - lambda||')

figure
A = load('output2_2.txt');
A = sort(A);
nevyazka = A(:, 1);
move = A(:, 2);
semilogy(move, nevyazka)
grid on
xlabel('delta')
ylabel('||Ax - lambda * x||')

figure
A = load('output2_3.txt');
A = sort(A);
x = A(:, 1);
move = A(:, 2);
plot(move, x)
grid on
xlabel('delta')
ylabel('||x* - x||')

figure
A = load('output3_1.txt');
A = sort(A);
lambda = A(:, 1);
move = A(:, 2);
plot(move, lambda)
grid on
xlabel('delta')
ylabel('||lambda* - lambda||')

figure
A = load('output3_2.txt');
A = sort(A);
nevyazka = A(:, 1);
move = A(:, 2);
plot(move, nevyazka)
grid on
xlabel('delta')
ylabel('||Ax - lambda * x||')

figure
A = load('output3_3.txt');
A = sort(A);
x = A(:, 1);
move = A(:, 2);
plot(move, x)
grid on
xlabel('delta')
ylabel('||x* - x||')