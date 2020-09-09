%% Demo on how to use LADEL (without permutation)
if exist('solver')
    solver.delete();
    clear solver
end

n = 10;
M = sprand(n,n, 1e-1, 1) + 2*speye(n);
M = (M+M')/2; %make sure M is symmetric
x = rand(n,1);

%% Example 1: factorize and solve
solver = ladel(n);
solver.factorize(M);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

%% Example 2: factorize_advanced and row_mod
Mbasis = sprand(n,n, 8e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;

% Make the n/2 row/column only contain a diagonal element
M(n/2,:) = zeros(1,n);
M(:,n/2) = zeros(n,1);
d = rand(1);
M(n/2,n/2) = 1; 

Mbasis = Mbasis + M; %make sure entries of M are in Mbasis

[L,D] = solver.factorize_advanced(M, Mbasis);
% LD = ldlchol(M);

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);
% y_chol = ldlsolve(LD, x);
% assert(norm(y_chol-M\x) < 1e-12);

%% ADD row using row_mod
[Lupd, Dupd] = solver.row_mod(n/2, Mbasis(:,n/2), full(Mbasis(n/2,n/2)));
% LD = ldlrowmod(LD, n/2, Mbasis(:,n/2));

Mupd = M;
Mupd(:,n/2) = Mbasis(:,n/2);
Mupd(n/2,:) = Mbasis(n/2,:);

% y_chol = ldlsolve(LD, x);
% assert(norm(y_chol-Mupd\x) < 1e-12);

y = solver.dense_solve(x);
assert(norm(y-Mupd\x) < 1e-12);

%% DELETE row using row_mod
% We delete the row we just added again

solver.row_mod(n/2);
y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

solver.delete();