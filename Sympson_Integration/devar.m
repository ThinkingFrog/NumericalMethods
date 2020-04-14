x = -1 : 0.0001 : 1;
y = 2520 * x.^2 - 480 * x - 240;
plot(x, y)
grid on
min = min(y)
max = max(y)