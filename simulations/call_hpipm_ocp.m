function [ x_hpipm, t, status, iter ] = call_hpipm_ocp(N, nx, nu, ng, A, B, Q, R, QT, F, f, xub, xlb, uub, ulb, tol, x0, x_warm)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    codegen_data = 1;

    dim = hpipm_ocp_qp_dim(N);
    dim.set('nx', nx, 0, N);
    dim.set('nu', nu, 0, N-1); % controls
    dim.set('nbx', nx, 0, N-1); % state bounds
    dim.set('nbu', nu, 0, N-1); % control bounds
    dim.set('ng', ng, N); % general linear constraints
    %dim.set('ns', ns, 0, N); % slacks

    % print to shell
    %dim.print_C_struct();
    % codegen
    if codegen_data
        dim.codegen('ocp_qp_data.c', 'w');
    end

    %%% qp %%%
    qp = hpipm_ocp_qp(dim);

    % dynamics
    qp.set('A', full(A), 0, N-1);
    qp.set('B', full(B), 0, N-1);

    % cost
    qp.set('Q', full(Q), 0, N-1);
    qp.set('S', zeros(nx, nu), 0, N-1);
    qp.set('R', full(R), 0, N-1);
    qp.set('Q', full(QT), N);
    qp.set('q', zeros(nx, 1), 0, N);
    qp.set('r', zeros(nu, 1), 0, N-1);

    % constraints
    Jbu = eye(nu);
    Jbx = eye(nx);
    qp.set('Jbx', Jbx, 0, N-1);

    qp.set('lbx', x0, 0);
    qp.set('ubx', x0, 0);
    qp.set('lbx', xlb, 1, N-1);
    qp.set('ubx', xub, 1, N-1);

    qp.set('Jbu', Jbu, 0, N-1);
    qp.set('lbu', ulb, 0, N-1);
    qp.set('ubu', uub, 0, N-1);

%         qp.set('lg', -1e20*ones(nf,1), N);
    qp.set('ug', f , N);
    qp.set('C', full(F), N);
    
    % codegen
    if codegen_data
        qp.codegen('ocp_qp_data.c', 'a');
    end

    %%% sol %%%
    sol = hpipm_ocp_qp_sol(dim);

    %%% solver arg %%%
    %mode = 'speed_abs';
%     mode = 'speed';
    mode = 'balance';
    %mode = 'robust';
    % create and set default arg based on mode
    arg = hpipm_ocp_qp_solver_arg(dim, mode);

    % overwrite default argument values
%     arg.set('mu0', 1e4);
%     arg.set('iter_max', 20);
%     arg.set('tol_stat', 1e-4);
%     arg.set('tol_eq', 1e-5);
%     arg.set('tol_ineq', 1e-5);
%     arg.set('tol_comp', 1e-5);
%     arg.set('reg_prim', 1e-12);
    %arg.set('warm_start', 1);

    % codegen
    if codegen_data
        arg.codegen('ocp_qp_data.c', 'a');
    end

    %%% solver %%%
    solver = hpipm_ocp_qp_solver(dim, arg);

    % arg which are allowed to be changed
    solver.set('iter_max', 10000);
    solver.set('tol_stat', tol);
    solver.set('tol_eq', tol);
    solver.set('tol_ineq', tol);
    solver.set('tol_comp', tol);

    % set solution guess
    if isempty(x_warm)
        sol.set('x', x0, 0);
        xk = x0;
        for k = 1:N
            xk = A*xk;
            sol.set('x', xk, k);
        end
        sol.set('u', zeros(nu,1), 0, N-1);
    else
        for k = 0:N-1
            sol.set('x', x_warm(k*(nu+nx)+1:k*nu+(k+1)*nx), k);
            sol.set('u', x_warm((k+1)*nx+k*nu+1:(k+1)*(nu+nx)), k);
%             x_hpipm(k*(nu+nx)+1:k*nu+(k+1)*nx) = x(k*nx+1:(k+1)*nx);
%             x_hpipm((k+1)*nx+k*nu+1:(k+1)*(nu+nx)) = u(k*nu+1:(k+1)*nu);
        end
        sol.set('x', x_warm(N*(nu+nx)+1:N*nu+(N+1)*nx), N);
    end
    
    % solve qp
    nrep = 1;
    tic
    for rep=1:nrep
        solver.solve(qp, sol);
    end
    solve_time = toc;

    % get solution statistics
    status = solver.get('status');
    t = solver.get('time_ext');
    iter = solver.get('iter');
%     res_stat = solver.get('max_res_stat');
%     res_eq = solver.get('max_res_eq');
%     res_ineq = solver.get('max_res_ineq');
%     res_comp = solver.get('max_res_comp');
%     stat = solver.get('stat');
%     if (~strcmp(travis_run, 'true'))
% %             status
% %             iter
% %             res_stat
% %             res_eq
% %             res_ineq
% %             res_comp
% %             fprintf('\nprint solver statistics\n');
%         fprintf('average solve time over %d runs: %e [s]\n', nrep, solve_time/nrep);
% %             fprintf('solve time of last run (measured in mex interface): %e [s]\n', time_ext);
% %             fprintf('iter\talpha_aff\tmu_aff\t\tsigma\t\talpha_prim\talpha_dual\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp\n');
% %             for ii=1:iter+1
% %                 fprintf('%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n', stat(ii,1), stat(ii,2), stat(ii,3), stat(ii,4), stat(ii,5), stat(ii,6), stat(ii,7), stat(ii,8), stat(ii,9), stat(ii,10), stat(ii,11));
% %             end
%     end



    % get solution
    % x
    x = sol.get('x', 0, N);
    x = reshape(x, nx, N+1);
    % u
    u = sol.get('u', 0, N-1);
    u = reshape(u, nu, N);
    % slack values
%         su = sol.get('su', 0, N);
%         sl = sol.get('sl', 0, N);
    x_hpipm = zeros(N*nu+(N+1)*nx,1);
    for k = 0:N-1
        x_hpipm(k*(nu+nx)+1:k*nu+(k+1)*nx) = x(k*nx+1:(k+1)*nx);
        x_hpipm((k+1)*nx+k*nu+1:(k+1)*(nu+nx)) = u(k*nu+1:(k+1)*nu);
    end
    x_hpipm(N*(nu+nx)+1:N*nu+(N+1)*nx) = x(N*nx+1:(N+1)*nx);

    
%     diff_sol = norm(x_qpalm - x_hpipm)/norm(x_qpalm)
% 
%     obj = 0.5*x_hpipm'*prob.Q*x_hpipm
%     max(max(prob.A*x_hpipm - ub, 0))
%     max(max(lb - prob.A*x_hpipm, 0))

end

