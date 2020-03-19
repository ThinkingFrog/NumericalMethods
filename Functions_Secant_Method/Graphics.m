p = [3, 4, -12, 0, 1];
x = -3 : 0.01 : 3;
y = polyval(p, x);
plot(x, y)
title('P(x) = 3x^4 + 4x^3 - 12x^2 + 1')
grid on
hold on
r = roots(p);
yr = polyval(p, r);
plot(r, yr, 'r*');

figure
y = f11(x);
plot(x, y);
title('F(x) = 6x + 3 - 5^x')
grid on
hold on
xr1 = fzero(@f11, -1);
yr1 = f11(xr1);
plot(xr1, yr1, 'r*');
xr2 = fzero(@f11, 1);
yr2 = f11(xr2);
plot(xr2, yr2, 'r*');