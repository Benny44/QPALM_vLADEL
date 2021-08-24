%Portfolio optimization
%Sparse regressor selection
clear; close all;

options.qpalm_matlab = false;
options.qpalm_c = true;
options.osqp = true;
options.qpoases = true;
options.gurobi = true;
options.VERBOSE = false;
options.SCALING_ITER = 10;
options.EPS_ABS = 1e-6;
options.TIME_LIMIT = 600;

Tqpalm_matlab = [];
Tqpalm_c = [];
Tosqp = [];
Tqpoases = [];
Tgurobi = [];
% rng(1)

nb_gamma = 10;
n_values = 100:20:1000;
nb_n = length(n_values);

options.return_solvers = true;

%Make use of special features of qpoases for this sequential problem
qpoases = options.qpoases;
options.qpoases = false;


for i = 1:nb_n
    n = n_values(i);
    m = n+1;
    
    k = ceil(n/10);
    F = sprandn(n,k,5e-1,1e-4);
    D = diag(rand(n,1)*sqrt(k));
    mu = randn(n,1);
%     mu(1:2:end) = mu(1:2:end)/1000;
    Q = blkdiag(D, eye(k));
    Q = sparse(Q);
    
%     m = m-1;
%     A = [F' -eye(k);
%         eye(n) zeros(n,k)];
%     A = sparse(A);
%     lb = [zeros(k+n,1)];
%     ub = [zeros(k,1); ones(n,1)];
    A = [ones(1,n) zeros(1,k); 
        F' -eye(k);
        eye(n) zeros(n,k)];
    A = sparse(A);
    lb = [1; zeros(k+n,1)];
    ub = [1; zeros(k,1); 1e20*ones(n,1)];

    
    prob.Q = Q; prob.A = A; prob.lb = lb; prob.ub = ub;
    
    qpalm_matlab_time = 0;
    qpalm_c_time = 0;
    osqp_time = 0;
    qpoases_time = 0;
    gurobi_time = 0;
    
    options.x = [];
    options.update = false;
    options.solvers_setup = false;
    options.return_solvers = true;
    
    gamma_range = logspace(-2,2,nb_gamma);
    
    for gamma = gamma_range;
        if gamma == max(gamma_range)
            options.return_solvers = false;
        end
        q = [-1/2/gamma*mu; zeros(k,1)];
        prob.q = q;
        [ X, timings, iter, status, options, stats ] = compare_QP_solvers(prob, options);
%         [X, timings, options] = compare_QP_solvers(prob, options);
        if options.qpalm_matlab , qpalm_matlab_time = qpalm_matlab_time + timings.qpalm_matlab; end
        if options.qpalm_c , qpalm_c_time = qpalm_c_time + timings.qpalm_c; end
        if options.osqp, osqp_time = osqp_time + timings.osqp; end
        if options.qpoases, qpoases_time = qpoases_time + timings.qpoases; end
        if options.gurobi, gurobi_time = gurobi_time + timings.gurobi; end
        
        %% Apply the first input (of the qpalm solution), shift the vectors for the next OCP and fix the bounds
        if strcmp(status.qpalm_c, 'solved')

            x_warm = X.qpalm_c;
            options.x = x_warm;
            options.update = true;
            options.solvers_setup = true;


        else
            error('Failed to solve for some reason');
        end
    end
    
    if qpoases
        qpoases_time = 0;
        q = [-1/2/gamma*mu; zeros(k,1)];
        qpoases_options = qpOASES_options('default', 'printLevel', 0, 'terminationTolerance', options.EPS_ABS, 'maxCpuTime', options.TIME_LIMIT, 'maxIter', 100000000);
        Iter_qpoases = [];

    %     [X_qpoases,fval,Status_qpoases,Iter_qpoases,lambda,auxOutput] = qpOASES(H,g',A,lb',ub',lbA',ubA',qpoases_options);
    %     [X_qpoases,fval,Status_qpoases,Iter_qpoases,lambda,auxOutput] = qpOASES(H,g',A,lb',ub',lbA',ubA',qpoases_options);
        gamma = gamma_range(1);
        [QP,x,fval,status_qpoases, iter_qpoases, ~, auxOutput] = qpOASES_sequence( 'i',Q,q,A,[],[],lb,ub,qpoases_options );
%         Status_qpoases{1} = status;
%         X_qpoases{1} = x;
%         Iter_qpoases = [Iter_qpoases iter];
        qpoases_time = qpoases_time + auxOutput.cpuTime;
        for gamma = gamma_range(2:end)
            q = [-1/2/gamma*mu; zeros(k,1)];
            [x,fval,status_qpoases, iter_qpoases, ~, auxOutput] = qpOASES_sequence( 'h',QP,q,[],[],lb,ub,qpoases_options );
%             Status_qpoases{i} = status;
%             X_qpoases{i} = x;
%             Iter_qpoases = [Iter_qpoases iter];
%             Tqpoases = [Tqpoases auxOutput.cpuTime];
            qpoases_time = qpoases_time + auxOutput.cpuTime;

        end
        qpOASES_sequence( 'c',QP );  
         
        Tqpoases(i) = qpoases_time/nb_gamma;
        
    end
    
    if options.qpalm_matlab, Tqpalm_matlab(i) = qpalm_matlab_time/nb_gamma; end
    if options.qpalm_c, Tqpalm_c(i) = qpalm_c_time/nb_gamma; end
    if options.osqp, Tosqp(i) = osqp_time/nb_gamma; end
    if options.qpoases, Tqpoases(i) = qpoases_time/nb_gamma; end
    if options.gurobi, Tgurobi(i) = gurobi_time/nb_gamma; end
    
    if options.qpalm_matlab, Iter_qpalm_matlab(i) = iter.qpalm_matlab; end
    if options.qpalm_c, Iter_qpalm_c(i) = iter.qpalm_c; end
    if options.osqp, Iter_osqp(i) = iter.osqp; end
    if options.qpoases, Iter_qpoases(i) = iter.qpoases; end
    if options.gurobi, Iter_gurobi(i) = iter.gurobi; end
    
    if options.qpalm_matlab, Status_qpalm_matlab{i} = status.qpalm_matlab; end
    if options.qpalm_c, Status_qpalm_c{i} = status.qpalm_c; end
    if options.osqp, Status_osqp{i} = status.osqp; end
    if options.qpoases, Status_qpoases{i} = status.qpoases; end
    if options.gurobi, Status_gurobi{i} = status.gurobi; end
    
    if options.qpalm_matlab, X_qpalm_matlab{i} = X.qpalm_matlab; end
    if options.qpalm_c, X_qpalm_c{i} = X.qpalm_c; end
    if options.osqp, X_osqp{i} = X.osqp; end
    if options.qpoases, X_qpoases{i} = X.qpoases; end
    if options.gurobi, X_gurobi{i} = X.gurobi; end
    
end

if options.osqp, options.osqp_solver.delete(); end
if options.qpalm_c, options.qpalm_solver.delete(); end

save('output/Portfolio_sequential')
% save('output/Portfolio', 'n_values','Tqpalm_matlab','Tqpalm_c','Tosqp','Tqpoases','Tgurobi');

%% Plot results

plot_portfolio_QP_comparison('output/Portfolio_sequential')
    
