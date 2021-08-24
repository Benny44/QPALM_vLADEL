function plot_MPC_QP_comparison( file, sequential )
%Helper function to plot results from simulation

load(file)
close all
figure

semilogy(n_values(1:length(Tqpalm_c)), Tqpalm_c, 'b',...
    n_values(1:length(Tosqp)), Tosqp, 'r',...
    n_values(1:length(Tqpoases)), Tqpoases, 'g',...
    n_values(1:length(Tgurobi)), Tgurobi, 'k',...
    n_values(1:length(Thpipm)), Thpipm, 'y');
%      
grid on
set(gca,'fontsize',14)
if sequential
    xlabel('Time step')
else
    xlabel('Number of primal variables')
end

ylabel('Runtime (s)')
% title(plot_title);
legend('QPALM', 'OSQP', 'qpOASES','Gurobi','HPIPM','Location','northwest')

end

