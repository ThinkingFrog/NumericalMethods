figure
A = load('output1_1.txt');
x = A(:, 1);
y = A(:, 2);
y = -log10(y);
plot(y, x)
grid on
ylabel('n')
xlabel('-log10(e)')

figure
A = load('output1_2.txt');
x = A(:, 1);
y = A(:, 2);
plot(x, y)
grid on
xlabel('x0')
ylabel('n')

figure
A = load('output1_3.txt');
x = A(:, 1);
y = A(:, 2);
y = -log10(y);
plot(y, x)
grid on
ylabel('n')
xlabel('-log10(e)')

figure
A = load('output1_4.txt');
x = A(:, 1);
y = A(:, 2);
plot(x, y)
grid on
xlabel('x0')
ylabel('n')

figure
A = load('output2_1.txt');
x = A(:, 1);
y = A(:, 2);
y = -log10(y);
plot(y, x)
grid on
ylabel('n')
xlabel('-log10(e)')

figure
A = load('output2_2.txt');
x = A(:, 1);
y = A(:, 2);
plot(x, y)
grid on
xlabel('x0')
ylabel('n')

figure
A = load('output2_3.txt');
x = A(:, 1);
y = A(:, 2);
y = -log10(y);
plot(y, x)
grid on
ylabel('n')
xlabel('-log10(e)')

figure
A = load('output2_4.txt');
x = A(:, 1);
y = A(:, 2);
plot(x, y)
grid on
xlabel('x0')
ylabel('n')