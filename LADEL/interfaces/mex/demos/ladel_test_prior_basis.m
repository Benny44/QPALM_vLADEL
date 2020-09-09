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

solver = ladel(n);

%% First use factorize_advanced with a basis matrix Mbasis
tic;
Mbasis = sprand(n,n, 2e-1, 1) + 3*speye(n);
Mbasis = (Mbasis+Mbasis')/2;

% Make sure the changed_row row/column only contain a diagonal element
M(changed_row,:) = zeros(nb_change,n);
M(:,changed_row) = zeros(n,nb_change);
M(changed_row,changed_row) = speye(nb_change); 
% Make sure entries of M are in Mbasis
Mbasis = Mbasis + M; 
Mbasis(changed_row, changed_row) = sparse(diag(rand(nb_change,1)));
times.generate.Mbasis = toc;

tic;
solver.factorize_advanced(M, Mbasis, ordering);
times.ladel.factorize = toc;

y = solver.dense_solve(x);
assert(norm(y-M\x) < TOL);


%% ADD rows using factorization_with_prior_basis
Mupd = M;
Mupd(:,changed_row) = Mbasis(:,changed_row);
Mupd(changed_row,:) = Mbasis(changed_row,:);

tic;
solver.factorize_with_prior_basis(Mupd);
times.ladel.factorize_with_prior_basis = toc;

y = solver.dense_solve(x);
assert(norm(y-Mupd\x) < TOL);


%reference time
tic;
solver.factorize_advanced(Mupd, Mbasis, ordering);
times.ladel.factorize_reference = toc;

solver.delete();

times.ladel
