%% Demo on how to use LADEL (with multiple rows)
if exist('solver')
    solver.delete();
    clear solver
end

test_cholmod = true;

ordering = 1; %1 for AMD, 0 for natural ordering

n = 1000;
nb_change = 100; %number of changing rows
TOL = 1e-10 * 2^(nb_change/2);
changed_row = randperm(n,nb_change);

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
assert(norm(y-M\x) < TOL);

%% Example 2: factorize_advanced and row_mod
tic;
Mbasis = sprand(n,n, 2e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;

% Make the changed_row row/column only contain a diagonal element
M(changed_row,:) = zeros(nb_change,n);
M(:,changed_row) = zeros(n,nb_change);
M(changed_row,changed_row) = speye(nb_change); 

Mbasis = Mbasis + M; %make sure entries of M are in Mbasis
Mbasis(changed_row, changed_row) = sparse(diag(rand(nb_change,1)));

times.generate.Mbasis = toc;

% clear all
% load('problematic');
% clear times
% clear solver
% solver = ladel(n);
% test_cholmod = true;

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
assert(norm(y-M\x) < TOL);
if test_cholmod
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-M\x) < TOL);
end

%% ADD row using row_mod
clear row;
row = Mbasis(:,changed_row);
d = diag(Mbasis);
d = full(d(changed_row));
tic;
solver.row_mod(changed_row, row, d);
times.ladel.rowadd = toc;

if test_cholmod
%     row = Mbasis(:,changed_row);
    tic;
    for k = changed_row
        row = Mbasis(LD_p, k);
        LD = ldlrowmod(LD, LD_pinv(k), row);
    end
    times.chol.rowadd = toc;
end

Mupd = M;
Mupd(:,changed_row) = Mbasis(:,changed_row);
Mupd(changed_row,:) = Mbasis(changed_row,:);

if test_cholmod
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-Mupd\x) < TOL);
end

y = solver.dense_solve(x);
assert(norm(y-Mupd\x) < TOL);

%% DELETE row using row_mod
% We delete the rows we just added again

if test_cholmod
    tic;
    for k = changed_row
        LD = ldlrowmod(LD, LD_pinv(k));
    end
    times.chol.rowdel = toc;
    y_chol = ldl_perm_solve(LD, LD_p, x);
    assert(norm(y_chol-M\x) < TOL);
end

tic;
solver.row_mod(changed_row);
times.ladel.rowdel = toc;

y = solver.dense_solve(x);
assert(norm(y-M\x) < TOL);

solver.delete();

times.ladel
if test_cholmod
    times.chol
end