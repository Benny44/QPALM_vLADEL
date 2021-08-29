function [results, options] = run_test_set( myFolder, options )
%Run test set in myFolder for solvers given in options
%Returns the options and a struct results which contains
%X (primal solutions), 
%T (timings), 
%I (iterations), 
%S (status),
%myFiles (list of filenames ran)
%Y (dual solutions), if options.return_y is set and only for qpalm and
%ipopt (needed for the cutest test set) 


if ~isfield(options, 'qpalm_c')
    options.qpalm_c = false;
end

if ~isfield(options, 'qpalm_matlab')
    options.qpalm_matlab = false;
end

if ~isfield(options, 'osqp')
    options.osqp = false;
end

if ~isfield(options, 'qpoases')
    options.qpoases = false;
end

if ~isfield(options, 'gurobi')
    options.gurobi = false;
end

if ~isfield(options, 'ipopt')
    options.ipopt = false;
end

if ~isfield(options, 'NewtonKKTqp')
    options.NewtonKKTqp = false;
end

if ~isfield(options, 'return_y')
    options.return_y = false;
end

current_path = fileparts(mfilename('fullpath'));
cd(current_path);

myFiles = dir(strcat(myFolder, '/*.mat'));
ll                   = length(myFiles);

% T.qpalm_matlab = [];
% T.qpalm_c = [];
% T.osqp = [];
% T.qpoases = [];
% T.gurobi = [];
% T.ipopt = [];
% 
% I.qpalm_matlab = [];
% I.qpalm_c = [];
% I.osqp = [];
% I.qpoases = [];
% I.gurobi = [];
% I.ipopt = [];
% 
% X.qpalm_matlab = [];
% X.qpalm_c = [];
% X.osqp = [];
% X.qpoases = [];
% X.gurobi = [];
% X.ipopt = [];
% 
% S.qpalm_matlab = [];
% S.qpalm_c = [];
% S.osqp = [];
% X.qpoases = [];
% X.gurobi = [];
% X.ipopt = [];

T = {};
I = {};
X = {};
S = {};


for i = 1:ll
    baseFileName = myFiles(i).name;
%     new{i}       = char(baseFileName(1:end-4));
    fullFileName = fullfile(baseFileName);
%     maros_files{i} = fullFileName;
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    matData{i}   = load(fullFileName); 
    if isfield(matData{i}, 'Data')  %Cutest format
        matData{i} = matData{i}.Data;
        P            = matData{i}.Q;
        l            = matData{i}.bl;
        u            = matData{i}.bu;
        lb           = matData{i}.cl;
        ub           = matData{i}.cu;
        prob.l = l; prob.u = u;
    else %Maros meszaros format
        n            = matData{i}.n;
        m            = matData{i}.m;
        P            = matData{i}.P;
        q            = matData{i}.q;
        lb           = matData{i}.l;
        ub           = matData{i}.u;
        
    end
    n            = matData{i}.n;
    m            = matData{i}.m;
    q            = matData{i}.q;
    A            = matData{i}.A;
    
    %delete empty rows of A (assuming lb <= 0 <= ub for these rows)
    b = A*rand(n,1);
    I = find(abs(b) == 0);
    A(I,:) = [];
    lb(I,:) = [];
    ub(I,:) = [];
    
    prob.Q = P; prob.q = q; prob.lb = lb; prob.ub = ub; prob.A = A;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
      
    [X_sol, timings, iter, status, options, ~] = compare_QP_solvers(prob, options);
    if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
    if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
    if options.osqp, osqp_time = osqp_time + timings.osqp; end
    if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
    if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
    
    if options.qpalm_matlab, T.qpalm_matlab(i) = qpalm_matlab_time; end
    if options.qpalm_c, T.qpalm_c(i) = qpalm_c_time; end
    if options.osqp, T.osqp(i) = osqp_time; end
    if options.qpoases, T.qpoases(i) = qpoases_time; end
    if options.gurobi, T.gurobi(i) = gurobi_time; end

    if options.qpalm_matlab, I.qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_c, I.qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, I.osqp(i) = iter.osqp; end
    if options.qpoases, I.qpoases(i) = iter.qpoases; end
    if options.gurobi, I.gurobi(i) = iter.gurobi; end
    
    if options.qpalm_matlab, S.qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, S.qpalm_c{i} = status.qpalm_c; end
    if options.osqp, S.osqp{i} = status.osqp; end
    if options.qpoases, S.qpoases{i} = status.qpoases; end
    if options.gurobi, S.gurobi{i} = status.gurobi; end
    
    if options.qpalm_matlab, X.qpalm_matlab{i} = X_sol.qpalm_matlab; end
    if options.qpalm_c, X.qpalm_c{i} = X_sol.qpalm_c; end
    if options.osqp, X.osqp{i} = X_sol.osqp; end
    if options.qpoases, X.qpoases{i} = X_sol.qpoases; end
    if options.gurobi, X.gurobi{i} = X_sol.gurobi; end
    
    if options.return_y
        if options.qpalm_c, Y.qpalm_c{i} = options.y.qpalm_c; end
        if options.ipopt, Y.ipopt{i} = options.y.ipopt; end
        options.y = [];
    end
end

results.X = X;
results.T = T;
results.I = I;
results.S = S;
results.myFiles = myFiles;
if options.return_y
    results.Y = Y;
end

end

