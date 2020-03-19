figure
A = load('Dependings1.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('eps')
ylabel('||x* - x||')

figure
A = load('Dependings2.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('eps')
ylabel('||Ax - b||')

figure
A = load('Dependings3.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('det(A)')
ylabel('||x* - x||')

figure
A = load('Dependings4.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('det(A)')
ylabel('||Ax - b||')

figure
A = load('Dependings5.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('eps')
ylabel('N')

figure
A = load('Dependings6.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('x0')
ylabel('N')

figure
A = load('Dependings7.txt');
A = sort(A);
d2 = A(:, 1);
d1 = A(:, 2);
loglog(d1, d2)
grid on
xlabel('det(A)')
ylabel('N')