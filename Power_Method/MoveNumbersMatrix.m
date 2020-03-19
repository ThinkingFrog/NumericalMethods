n = 12;
for i = 0 : 80
x=[(1 + i) 8 16 24 32 40 48 56 64 72 80 88];
A=diag(x);
w=[1 1 1 1 1 1 1 1 1 1 1 1];
w=w/norm(w);
w=w';
E=eye(n);
P=E-2/(w'*w)*(w*w');
Q=inv(P);
newA=Q*A*P;
fname = sprintf ('%s %i', 'Matrix', i);
fname = sprintf ('%s %s', fname, '.txt');
dlmwrite(fname, newA, 'precision', 16, 'delimiter', ' ')
cond(A);
y=det(newA);

A = load('Matrix.txt');
[R, D] = eig(A);
R = R';

lambda = max(max(D));

end