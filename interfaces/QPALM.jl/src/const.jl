const QPALM_INFTY = 1e20

const SOLUTION_PRESENT = [:Solved, :Max_iter_reached]

const status_map = Dict{Int,Symbol}(
     0  => :Error,
     1  => :Solved,
     2  => :Dual_terminated,
    -2  => :Max_iter_reached,
    -3  => :Primal_infeasible,
    -4  => :Dual_infeasible,
    -5  => :Time_limit_reached,
    -10 => :Unsolved
)
