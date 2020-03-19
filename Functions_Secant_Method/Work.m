x = [-3 -2.5];
eps=1e-15;
options = optimset('Display','iter','TolX', eps);
x1=fzero(@f12,x,options);
x = [1.5 2];
x2=fzero(@f11,x,options);
dlmwrite('result1.txt', x1, 'precision', 16)
dlmwrite('result2.txt', x2, 'precision', 16)