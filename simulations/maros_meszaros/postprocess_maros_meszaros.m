function [] = postprocess_maros_meszaros(results, options)

close all

current_path = fileparts(mfilename('fullpath'));
cd(current_path);

TIME_LIMIT = options.TIME_LIMIT;

Tosqp = results.T.osqp;
Tqpalm_c = results.T.qpalm_c;
Tgurobi = results.T.gurobi;

Status_osqp = results.S.osqp;
Status_qpalm_c = results.S.qpalm_c;
Status_gurobi = results.S.gurobi;

% load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/MM_OSQP.mat')
[gs_osqp, fail_rate_osqp, Tosqp] = compute_geometric_mean(Tosqp, Status_osqp, 'solved', TIME_LIMIT);

% load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/MM_GUROBI.mat')
[gs_gurobi, fail_rate_gurobi, Tgurobi] = compute_geometric_mean(Tgurobi, Status_gurobi, 'OPTIMAL', TIME_LIMIT);

% load('/home/ben/Documents/Projects/QPALM/simulations/results/journal_paper/MM_QPALM_C.mat')
[gs_qpalm, fail_rate_qpalm, Tqpalm_c] = compute_geometric_mean(Tqpalm_c, Status_qpalm_c, 'solved', TIME_LIMIT);

cd('../results/journal_paper/table_MM/');
fid = fopen('MM.tex', 'w');
fprintf(fid, 'Runtime (sgm) & %.4f & %.4f & %.4f\\\\\n', gs_qpalm, gs_osqp, gs_gurobi);
fprintf(fid, 'Failure rate [\\%%] & %.4f & %.4f & %.4f\n', fail_rate_qpalm, fail_rate_osqp, fail_rate_gurobi);
fcl = fclose(fid);


if options.EPS_ABS == 1e-3, mytol = '1e3'; end
if options.EPS_ABS == 1e-6, mytol = '1e6'; end
mytexfile = strcat(' table_MM_', mytol);
mytexcmd = strcat('pdflatex', mytexfile);
cmd = system(mytexcmd);

cd(current_path);

%% Performance profiles

Tosqp(Tosqp==TIME_LIMIT) = inf;
Tqpalm_c(Tqpalm_c==TIME_LIMIT) = inf;
Tgurobi(Tgurobi==TIME_LIMIT) = inf;

%Compare QPALM and Gurobi
Tmin = min(Tqpalm_c, Tgurobi);
r_qpalm_c = Tqpalm_c./Tmin;
r_gurobi = Tgurobi./Tmin;

rmax = max(max(r_gurobi(r_gurobi~=inf & ~isnan(r_gurobi))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
% if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
% end


[xgu, ygu, fgu] = make_performance_profile(r_gurobi);
% if fgu
    xgu = [xgu rmax];
    ygu = [ygu ygu(end)];
% end

figure
plot(log10(xqp), yqp, 'b', log10(xgu) ,ygu, 'k')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'Gurobi','Location', 'SouthEast')

%Compare QPALM and OSQP
Tmin = min(Tqpalm_c, Tosqp);
r_qpalm_c = Tqpalm_c./Tmin;
r_osqp = Tosqp./Tmin;

rmax = max(max(r_osqp(r_osqp~=inf & ~isnan(r_osqp))), max(r_qpalm_c(r_qpalm_c~=inf & ~isnan(r_qpalm_c)))); 

[xqp, yqp, fqp] = make_performance_profile(r_qpalm_c);
% if fqp
    xqp = [xqp rmax];
    yqp = [yqp yqp(end)];
% end


[xos, yos, fos] = make_performance_profile(r_osqp);
% if fos
    xos = [xos rmax];
    yos = [yos yos(end)];
% end

figure
plot(log10(xqp), yqp, 'b', log10(xos) ,yos, 'r')
set(gca,'fontsize',14)
xlabel('log_{10}(f)')
ylabel('fraction of solver within f of best')
legend('QPALM', 'OSQP','Location', 'SouthEast')

end