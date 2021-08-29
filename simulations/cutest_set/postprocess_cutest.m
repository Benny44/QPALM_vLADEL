function [] = postprocess_cutest(results, options)

current_path = fileparts(mfilename('fullpath'));
cd(current_path);

close all

TIME_LIMIT = options.TIME_LIMIT;
Tqpalm_c = results.T.qpalm_c;
Tipopt = results.T.ipopt;

Status_qpalm_c = results.S.qpalm_c;
Status_ipopt = results.S.ipopt;

X_qpalm_c = results.X.qpalm_c;
X_ipopt = results.X.ipopt;

Y_qpalm_c = results.Y.qpalm_c;
Y_ipopt = results.Y.ipopt;

myFiles = results.myFiles;
for i = 1:length(myFiles)
    myFilenames{i} = myFiles(i).name;
end

%Infeasible problems
I = [];
index = find(strcmp(myFilenames, 'A2NNDNIL.mat')); I = [I index]; %PI
index = find(strcmp(myFilenames, 'A5NNDNIL.mat')); I = [I index]; %PI
index = find(strcmp(myFilenames, 'LINCONT.mat')); I = [I index]; %PI
index = find(strcmp(myFilenames, 'NASH.mat')); I = [I index]; %PI
index = find(strcmp(myFilenames, 'STATIC3.mat')); I = [I index]; %DI

Tqpalm_c(I) = [];
Status_qpalm_c(I) = [];
X_qpalm_c(I) = [];
Y_qpalm_c(I) = [];

Tipopt(I) = [];
Status_ipopt(I) = [];
X_ipopt(I) = [];
Y_ipopt(I) = [];

files = myFilenames;
files(I) = [];


%% Geometric average runtime and fail rates
[gs_qpalm, fail_rate_qpalm, Tqpalm_c] = compute_geometric_mean(Tqpalm_c, Status_qpalm_c, 'solved', TIME_LIMIT);
fail_rate_qpalm = fail_rate_qpalm*length(Tqpalm_c)/(length(Tqpalm_c) + 5); %correct for infeasible problems

[gs_ipopt, fail_rate_ipopt, Tipopt] = compute_geometric_mean(Tipopt, Status_ipopt, 'Solve_Succeeded', TIME_LIMIT);
fail_rate_ipopt = fail_rate_ipopt*length(Tipopt)/(length(Tipopt) + 5); %correct for infeasible problems
 

%% Performance profile and objective comparisons
Tqpalm_c(~strcmp(Status_qpalm_c, 'solved')) = inf; 
Tipopt(~strcmp(Status_ipopt, 'Solve_Succeeded')) = inf;

Obj_qpalm = [];
Obj_ipopt = [];
Tqp = [];
Tip = [];

x_tol = 1e-6;
nb_same = 0;
nb_dif = 0;
nb_both_fail = 0;
nb_qpalm_wins = 0;
nb_ipopt_wins = 0;

for i = 1:length(Tqpalm_c)
    xqp = X_qpalm_c{i};
    xip = full(X_ipopt{i});
    load(files{i});
    obj_qp = 0.5*xqp'*Data.Q*xqp + Data.q'*xqp;
    obj_ip = 0.5*xip'*Data.Q*xip + Data.q'*xip;
    feas_qp = norm([max(Data.A*xqp - Data.cu, 0)+min(Data.A*xqp - Data.cl, 0); max(xqp - Data.bu, 0)+min(xqp - Data.bl, 0)]);
    feas_ip = norm([max(Data.A*xip - Data.cu, 0)+min(Data.A*xip - Data.cl, 0); max(xip - Data.bu, 0)+min(xip - Data.bl, 0)]);
    diff_x = norm(xip-xqp)/norm(xqp);    
    if (Tqpalm_c(i) == inf && Tipopt(i) == inf)
        if (feas_qp < feas_ip && obj_qp < obj_ip)
            nb_qpalm_wins = nb_qpalm_wins + 1;
        elseif (feas_qp < feas_ip && obj_qp < obj_ip)
            nb_ipopt_wins = nb_ipopt_wins + 1;
        else
%             nb_same = nb_same + 1;
%             Tqp(nb_same) = Tqpalm_c(i);
%             Tip(nb_same) = Tipopt(i);
        end
    elseif diff_x < x_tol
        nb_same = nb_same + 1;
        Tqp(nb_same) = Tqpalm_c(i);
        Tip(nb_same) = Tipopt(i);
    else
        nb_dif = nb_dif + 1;
        Obj_qpalm(nb_dif) = obj_qp;
        if (Tqpalm_c(i) ~= inf && Tipopt(i) == inf)
            Obj_ipopt(nb_dif) = inf;
        else
            Obj_ipopt(nb_dif) = obj_ip;
        end
    end
end

%Compare objectives for different solutions (don't count the failures to solve)
nb_qpalm_wins = sum(Obj_qpalm < Obj_ipopt);
nb_ipopt_wins = sum(Obj_qpalm > Obj_ipopt);

%Compare timings for same solution of QPALM and Ipopt
Tmin = min(Tqp, Tip);
r_qpalm_c = Tqp./Tmin;
r_ipopt = Tip./Tmin;

rmax = max(max(r_ipopt(r_ipopt~=inf & ~isnan(r_ipopt))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xgu, ygu, fgu] = make_performance_profile(r_ipopt);
% if fgu
    xgu = [xgu rmax];
    ygu = [ygu ygu(end)];
% end


[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
% if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
% end

figure
plot(log10(xqp), yqp, 'b', log10(xgu) ,ygu, 'k')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'IPOPT','Location', 'SouthEast')

Tip(Tip==inf) = TIME_LIMIT;
Tqp(Tqp==inf) = TIME_LIMIT;
n = length(Tip);
all_success = {};
for k = 1:n
    all_success{k} = 'success';
end
[gs_ipopt, ~, ~] = compute_geometric_mean(Tip, all_success, 'success', TIME_LIMIT);

[gs_qpalm, ~, ~] = compute_geometric_mean(Tqp, all_success, 'success', TIME_LIMIT);


%% Compare second-order necessary conditions violations
nb_qpalm_dead = 0;
nb_ipopt_dead = 0;

for i = 1:length(Tqpalm_c)
    xqp = X_qpalm_c{i};
    xip = full(X_ipopt{i});
    yqp = Y_qpalm_c{i};
    yip = full(Y_ipopt{i});
    load(files{i});
    A = Data.A;
    n = Data.n;
    lb = Data.cl;
    ub = Data.cu;

    %Delete empty rows in A
    b = A*rand(n,1);
    I = find(abs(b) == 0);
    A(I,:) = [];
    lb(I,:) = [];
    ub(I,:) = [];

    A = [A; speye(Data.n)];

%     A = [A; -A];
    bu = [ub; Data.bu];
    bl = [lb; Data.bl];
    

    if Tqpalm_c(i) ~= inf && check_dead_point(xqp, yqp, A, bu, bl, Data.Q, options.EPS_ABS)
        nb_qpalm_dead = nb_qpalm_dead + 1;
        fprintf('QPALM dead at i = %d\n', i);
    end
    
    if Tipopt(i) ~= inf && check_dead_point(xip, yip, A, bu, bl, Data.Q, options.EPS_ABS)
        nb_ipopt_dead = nb_ipopt_dead + 1;
        fprintf('IPOPT dead at i = %d\n', i);
    end   
    
end


%% Summary table

fid = fopen('../results/journal_paper/table_cutest_1/cutest_small.tex', 'w');
fprintf(fid, 'Runtime (sgm) & %.4f & %.4f\\\\\n', gs_qpalm, gs_ipopt);
fprintf(fid, 'Optimal & %4d & %4d\\\\\n', nb_qpalm_wins, nb_ipopt_wins);
fprintf(fid, 'Dead points & %4d & %4d \\\\\n', nb_qpalm_dead, nb_ipopt_dead);
fprintf(fid, 'Failure rate [\\%%] & %.4f & %.4f \n', fail_rate_qpalm, fail_rate_ipopt);

fcl = fclose(fid);
cd('../results/journal_paper/table_cutest_1/');
mytexfile = ' table_cutest';
mytexcmd = strcat('pdflatex', mytexfile);
cmd = system(mytexcmd);

cd(current_path)

%% Big table with timings per problem

Tqpalm_c = results.T.qpalm_c;
Tipopt = results.T.ipopt;
Status_qpalm_c = results.S.qpalm_c;
Status_ipopt = results.S.ipopt;
X_qpalm_c = results.X.qpalm_c;
X_ipopt = results.X.ipopt;

myFiles = results.myFiles;
for i = 1:length(myFiles)
    files{i} = myFiles(i).name;
end

fid = fopen('../results/journal_paper/table_cutest_2/cutest.tex', 'w');

i = find(Tipopt > TIME_LIMIT);
if (~isempty(i))
    for j = i
        Status_ipopt{j} = 'Maximum_CpuTime_Exceeded';
    end
end

[~, nb_files] = size(files);

for i = 1:nb_files
    load(files{i});
    
    fprintf(fid, '%s & %5ld & %5ld & ', Data.name, Data.n, Data.m);
    
    if (strcmp(Status_qpalm_c{i}, 'solved'))
        qpalm_status = 'S';
    elseif (strcmp(Status_qpalm_c{i}, 'primal infeasible'))
        qpalm_status = 'PI';
    elseif (strcmp(Status_qpalm_c{i}, 'dual infeasible'))
        qpalm_status = 'DI';
    else
        qpalm_status = 'F';
    end
    
%     fprintf(fid, '%2s & ', qpalm_status);
    
    if (~strcmp(qpalm_status, 'S'))
        qpalm_time = qpalm_status;
        qpalm_obj = '/';
        fprintf(fid, '%s & ', qpalm_time);
    else
        qpalm_time = Tqpalm_c(i);
%         try
        qpalm_obj = 0.5*X_qpalm_c{i}'*Data.Q*X_qpalm_c{i} + Data.q'*X_qpalm_c{i};
%         catch
%             set_break = 1;
%         end
        fprintf(fid, '%.2e & ', qpalm_time);
    end
    
    if (strcmp(Status_ipopt{i}, 'Solve_Succeeded'))
        ipopt_status = 'S';
    elseif (strcmp(Status_ipopt{i}, 'Infeasible_Problem_Detected'))
        ipopt_status = 'PI';
    elseif (strcmp(Status_ipopt{i}, 'Diverging_Iterates'))
        ipopt_status = 'DI';
    else
        ipopt_status = 'F';
    end
    
%     fprintf(fid, '%2s & ', ipopt_status);
    
    if (~strcmp(ipopt_status, 'S'))
        ipopt_time = ipopt_status;
        ipopt_obj = '/';
        fprintf(fid, '%s & ', ipopt_time);
    else
        ipopt_time = Tipopt(i);
        ipopt_obj = full(0.5*X_ipopt{i}'*Data.Q*X_ipopt{i} + Data.q'*X_ipopt{i});
        fprintf(fid, '%.2e & ', ipopt_time);
    end
    
    if (strcmp(qpalm_status, 'S'))
        if (qpalm_obj >= 0)
            fprintf(fid, '\\phantom{-}');
        end
        fprintf(fid, '%.2e & ', qpalm_obj);
    else
        fprintf(fid, '%s & ', qpalm_obj);
    end
    
    if (strcmp(ipopt_status, 'S'))
        if (ipopt_obj >= 0)
            fprintf(fid, '\\phantom{-}');
        end
        fprintf(fid, '%.2e', ipopt_obj);
    else
        fprintf(fid, '%s', ipopt_obj);
    end
    
    fprintf(fid, '\\\\ \n');     
end

fcl = fclose(fid);
cd('../results/journal_paper/table_cutest_2/');
mytexfile = ' table_cutest';
mytexcmd = strcat('pdflatex', mytexfile);
cmd = system(mytexcmd);

cd(current_path);


