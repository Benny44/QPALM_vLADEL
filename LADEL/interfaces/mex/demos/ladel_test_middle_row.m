%% Demo on how to use LADEL (with permutation)
if exist('solver')
    solver.delete();
    clear solver
end

test_cholmod = true;

ordering = 0; %1 for AMD, 0 for natural ordering

n = 500;
changed_row = n/2;

tic
M = sprand(n,n, 1e-1, 1) + 2*speye(n);
M = (M+M')/2;
x = rand(n,1);
times.generate.M = toc;

%% Example 1: factorize and solve
solver = ladel(n);
% [L,D,p] = solver.factorize(M, ordering);
solver.factorize(M, ordering);

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

%% Example 2: factorize_advanced and row_mod
tic;
Mbasis = sprand(n,n, 2e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;

% Make the changed_row row/column only contain a diagonal element
M(changed_row,:) = zeros(1,n);
M(:,changed_row) = zeros(n,1);
d = rand(1);
M(changed_row,changed_row) = 1; 

Mbasis = Mbasis + M; %make sure entries of M are in Mbasis
times.generate.Mbasis = toc;

% clear all
% load('problematic');
% clear times
% clear solver
% solver = ladel(n);
% ordering = 1;
tic;
% [L,D,p] = solver.factorize_advanced(M, Mbasis, ordering);
solver.factorize_advanced(M, Mbasis, ordering);
times.ladel.factorize = toc;
if test_cholmod
    tic;
    [LD,~, LD_p] = ldlchol(M);
    LD_pinv = 1:n;
    LD_pinv(LD_p) = LD_pinv;
    times.chol.factorize = toc;
end

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);
if test_cholmod
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-M\x) < 1e-12);
end

%% ADD row using row_mod
clear row;
row = Mbasis(:,changed_row);
tic;
% [Lupd,Dupd,p] = solver.row_mod(changed_row, row, full(Mbasis(changed_row,changed_row)));
solver.row_mod(changed_row, row, full(Mbasis(changed_row,changed_row)));
times.ladel.rowadd = toc;

if test_cholmod
    row = Mbasis(:,changed_row);
    tic;
    row = row(LD_p);
    LD = ldlrowmod(LD, LD_pinv(changed_row), row);
    times.chol.rowadd = toc;
end

Mupd = M;
Mupd(:,changed_row) = Mbasis(:,changed_row);
Mupd(changed_row,:) = Mbasis(changed_row,:);

if test_cholmod
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-Mupd\x) < 1e-12);
end

y = solver.dense_solve(x);
assert(norm(y-Mupd\x) < 1e-12);

%% DELETE row using row_mod
% We delete the row we just added again

if test_cholmod
    tic;
    LD = ldlrowmod(LD, LD_pinv(changed_row));
    times.chol.rowdel = toc;
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-M\x) < 1e-12);
end

tic;
% [Lfinal, Dfinal, p] = solver.row_mod(changed_row);
solver.row_mod(changed_row);
times.ladel.rowdel = toc;

y = solver.dense_solve(x);
assert(norm(y-M\x) < 1e-12);

solver.delete();

times.ladel
times.chol