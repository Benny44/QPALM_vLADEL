function plot_portfolio_QP_comparison( file )
%Helper function to plot results from simulation

load(file)
close all
figure

semilogy(n_values(1:length(Tqpalm_c)), Tqpalm_c, 'b',...
    n_values(1:length(Tosqp)), Tosqp, 'r',...
    n_values(1:length(Tqpoases)), Tqpoases, 'g',...
    n_values(1:length(Tgurobi)), Tgurobi, 'k');
%      
grid on
set(gca,'fontsize',14)
% xlabel('Number of primal variables')
xlabel('Time step')
ylabel('Runtime (s)')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','Location','northwest')

end

