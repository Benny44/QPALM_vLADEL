using SparseArrays
linsys = "ladel" #"ladel" or "cholmod"

if linsys=="cholmod"
    using SuiteSparse: CHOLMOD
end

const Cc_int = if Sys.WORD_SIZE == 64
    Clonglong
else
    Cint
end

const Cc_uint = if Sys.WORD_SIZE == 64
    Culong
else
    Cuint
end


struct Carray_element
    x::Cdouble
    i::Cc_uint
end

struct ladel_sparse_matrix
    nzmax::Cc_int
    nrow::Cc_int
    ncol::Cc_int
    p::Ptr{Cc_int}
    i::Ptr{Cc_int}
    x::Ptr{Cdouble}
    nz::Ptr{Cc_int}
    values::Cc_int
    symmetry::Cc_int
end

if linsys=="ladel"
    sparse_matrix_ptr = Ref{QPALM.ladel_sparse_matrix}
elseif linsys=="cholmod"
    sparse_matrix_ptr = Ptr{CHOLMOD.C_Sparse}
else
    error("Unrecognised linsys")
end

struct Data
    n::Cc_uint
    m::Cc_uint
    Q::sparse_matrix_ptr
    A::sparse_matrix_ptr
    q::Ptr{Cdouble}
    c::Cdouble
    bmin::Ptr{Cdouble}
    bmax::Ptr{Cdouble}
end

struct Solution
    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
end

struct CInfo
    iter::Cc_int
    iter_out::Cc_int
    status::NTuple{32,Cchar}
    status_val::Cc_int
    pri_res_norm::Cdouble
    dua_res_norm::Cdouble
    dua2_res_norm::Cdouble
    objective::Cdouble
    dual_objective::Cdouble
    setup_time::Cdouble
    solve_time::Cdouble
    run_time::Cdouble
end

struct Settings
    max_iter::Cc_int
    inner_max_iter::Cc_int
    eps_abs::Cdouble
    eps_rel::Cdouble
    eps_abs_in::Cdouble
    eps_rel_in::Cdouble
    rho::Cdouble
    eps_prim_inf::Cdouble
    eps_dual_inf::Cdouble
    theta::Cdouble
    delta::Cdouble
    sigma_max::Cdouble
    sigma_init::Cdouble
    proximal::Cc_int
    gamma_init::Cdouble
    gamma_upd::Cdouble
    gamma_max::Cdouble
    scaling::Cc_int
    nonconvex::Cc_int
    verbose::Cc_int
    print_iter::Cc_int
    warm_start::Cc_int
    reset_newton_iter::Cc_int
    enable_dual_termination::Cc_int
    dual_objective_limit::Cdouble
    time_limit::Cdouble
    ordering::Cc_int
    factorization_method::Cc_int
    max_rank_update::Cc_int
    max_rank_update_fraction::Cc_int
end

struct Workspace
    data::Ptr{QPALM.Data}

    x::Ptr{Cdouble}
    y::Ptr{Cdouble}
    Ax::Ptr{Cdouble}
    Qx::Ptr{Cdouble}
    Aty::Ptr{Cdouble}
    x_prev::Ptr{Cdouble}
    initialized::Cc_int

    temp_m::Ptr{Cdouble}
    temp_n::Ptr{Cdouble}
    sigma::Ptr{Cdouble}
    sigma_inv::Ptr{Cdouble}
    sqrt_sigma_max::Cdouble
    nb_sigma_changed::Cc_int
    gamma::Cdouble
    gamma_maxed::Cc_int
    Axys::Ptr{Cdouble}
    z::Ptr{Cdouble}
    pri_res::Ptr{Cdouble}
    pri_res_in::Ptr{Cdouble}
    yh::Ptr{Cdouble}
    Atyh::Ptr{Cdouble}
    df::Ptr{Cdouble}
    x0::Ptr{Cdouble}
    xx0::Ptr{Cdouble}
    dphi::Ptr{Cdouble}
    neg_dphi::Ptr{Cdouble}
    dphi_prev::Ptr{Cdouble}
    d::Ptr{Cdouble}

    tau::Cdouble
    Qd::Ptr{Cdouble}
    Ad::Ptr{Cdouble}
    sqrt_sigma::Ptr{Cdouble}
    sqrt_delta::Cdouble
    eta::Cdouble
    beta::Cdouble
    delta::Ptr{Cdouble}
    alpha::Ptr{Cdouble}
    temp_2m::Ptr{Cdouble}
    delta2::Ptr{Cdouble}
    delta_alpha::Ptr{Cdouble}
    s::Ptr{Carray_element}
    index_L::Cc_int
    index_P::Cc_int
    index_J::Cc_int

    eps_pri::Cdouble
    eps_dua::Cdouble
    eps_dua_in::Cdouble
    eps_abs_in::Cdouble
    eps_rel_in::Cdouble

    delta_y::Ptr{Cdouble}
    Atdelta_y::Ptr{Cdouble}

    delta_x::Ptr{Cdouble}
    Qdelta_x::Ptr{Cdouble}
    Adelta_x::Ptr{Cdouble}

    D_temp::Ptr{Cdouble}
    E_temp::Ptr{Cdouble}

    solver::Ptr{Nothing}
    settings::Ptr{QPALM.Settings}
    scaling::Ptr{Nothing}
    solution::Ptr{QPALM.Solution}
    info::Ptr{QPALM.CInfo}

    timer::Ptr{Nothing}
end
