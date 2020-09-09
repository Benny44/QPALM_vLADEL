import Base.GC: @preserve

mutable struct Model
    workspace::Ptr{QPALM.Workspace}

    function Model()
        model = new(C_NULL)
        finalizer(QPALM.cleanup!, model)
        return model
    end
end

function Settings(model::Maybe{QPALM.Model} = nothing; kwargs...)
    for k in keys(kwargs)
        @assert k in fieldnames(QPALM.Settings) "Unknown setting: $k"
    end

    default_settings =
        if model == nothing
            ds = Ref{QPALM.Settings}()
            ccall(
                (:qpalm_set_default_settings, LIBQPALM_PATH),
                Nothing,
                (Ref{QPALM.Settings},),
                ds
            )
            ds[]
        else
            workspace = unsafe_load(model.workspace)
            unsafe_load(workspace.settings)
        end

    merged_settings = [
        k in keys(kwargs) ? kwargs[k] : getfield(default_settings, k)
        for k in fieldnames(QPALM.Settings)
    ]

    return QPALM.Settings(merged_settings...)
end

function to_sparse!(
    M::Maybe{AbstractMatrix} = nothing
)
    if linsys=="ladel"
        Sparse_M = QPALM.ladel_sparse_matrix(
            length(M.rowval),
            M.m,
            M.n,
            pointer(M.colptr),
            pointer(M.rowval),
            pointer(M.nzval),
            C_NULL,
            true,
            0
        )
    elseif linsys=="cholmod"
        Sparse_M = CHOLMOD.Sparse(M)
    end
    return Sparse_M
end


function setup!(
    model::QPALM.Model;
    Q::Maybe{AbstractMatrix} = nothing,
    q::Maybe{Vector{Float64}} = nothing,
    A::Maybe{AbstractMatrix} = nothing,
    bmin::Maybe{Vector{Float64}} = nothing,
    bmax::Maybe{Vector{Float64}} = nothing,
    kwargs...
)
    # Check problem dimensions
    if Q == nothing
        if q != nothing
            n = length(q)
        elseif A != nothing
            n = size(A, 2)
        else
            error("The problem does not have any variables!")
        end

    else
        n = size(Q, 1)
    end

    if A == nothing
        m = 0
    else
        m = size(A, 1)
    end


    # Check if parameters are nothing
    if ((A == nothing) & ( (bmin != nothing) | (bmax != nothing))) |
        ((A != nothing) & ((bmin == nothing) & (bmax == nothing)))
        error("A must be supplied together with bmin and bmax")
    end

    if (A != nothing) & (bmin == nothing)
        bmin = -Inf * ones(m)
    end
    if (A != nothing) & (bmax == nothing)
        bmax = Inf * ones(m)
    end

    if Q == nothing
        Q = sparse([], [], Float64[], n, n)
    end
    if q == nothing
        q = zeros(n)
    end
    if A == nothing
        A = sparse([], [], Float64[], m, n)
        bmin = zeros(m)
        bmax = zeros(m)
    end


    # Check if dimensions are correct
    if length(q) != n
        error("Incorrect dimension of q")
    end
    if length(bmin) != m
        error("Incorrect dimensions of bmin")
    end
    if length(bmax) != m
        error("Incorrect dimensions of bmax")
    end


    # Check or sparsify matrices
    if !issparse(Q)
        @warn("Q is not sparse. Sparsifying it now (it might take a while)")
        Q = sparse(Q)
    end
    if !issparse(A)
        @warn("A is not sparse. Sparsifying it now (it might take a while)")
        A = sparse(A)
    end

    # Convert lower and upper bounds from Julia infinity to QPALM infinity
    bmin = max.(bmin, -QPALM_INFTY)
    bmax = min.(bmax, QPALM_INFTY)

    Sparse_Q = to_sparse!(Q)
    Sparse_A = to_sparse!(A)
    
    settings = QPALM.Settings(; kwargs...)

    if linsys=="ladel"
        data = ccall(
            (:qpalm_julia_set_data, LIBQPALM_PATH),
            Ref{QPALM.Data},
            (Cc_int, 
             Cc_int, 
             Ref{QPALM.ladel_sparse_matrix},
             Ref{QPALM.ladel_sparse_matrix},
             Ptr{Cfloat},
             Cfloat,
             Ptr{Cfloat},
             Ptr{Cfloat}),
             n, m, Ref(Sparse_Q), Ref(Sparse_A), 
             pointer(q), 0, pointer(bmin), pointer(bmax)
             )

        model.workspace = ccall(
            (:qpalm_setup, LIBQPALM_PATH),
            Ptr{QPALM.Workspace},
            (Ref{QPALM.Data}, Ptr{QPALM.Settings}),
            data, Ref(settings)
            )
    elseif linsys=="cholmod"
        data = QPALM.Data(
            n, m,
            pointer(Sparse_Q),
            pointer(Sparse_A),
            pointer(q),
            0,
            pointer(bmin), pointer(bmax)
        )
        model.workspace = ccall(
            (:qpalm_setup, LIBQPALM_PATH),
            Ptr{QPALM.Workspace},
            (Ptr{QPALM.Data}, Ptr{QPALM.Settings}),
            Ref(data), Ref(settings)
        )
    end

    if model.workspace == C_NULL
        error("Error in QPALM setup")
    end
end

mutable struct Info
    iter::Int64
    iter_out::Int64
    status::Symbol
    status_val::Int64
    pri_res_norm::Float64
    dua_res_norm::Float64
    dua2_res_norm::Float64
    objective::Float64
    dual_objective::Float64
    setup_time::Float64
    solve_time::Float64
    run_time::Float64

    Info() = new()
end

function update!(
    model::QPALM.Model;
    q::Maybe{Vector{Float64}} = nothing,
    bmin::Maybe{Vector{Float64}} = nothing,
    bmax::Maybe{Vector{Float64}} = nothing,
    kwargs...
)
    if q != nothing
        ccall(
            (:qpalm_update_q, LIBQPALM_PATH),
            Cvoid,
            (Ptr{QPALM.Workspace}, Ptr{Cfloat}),
            model.workspace, pointer(q),
        )
    end

    # Convert lower and upper bounds from Julia infinity to QPALM infinity
    if bmin != nothing
        bmin = max.(bmin, -QPALM_INFTY)
    end
    if bmax != nothing
        bmax = min.(bmax, QPALM_INFTY)
    end
    
    if bmin != nothing || bmax != nothing
        ccall(
            (:qpalm_update_bounds, LIBQPALM_PATH),
            Cvoid,
            (Ptr{QPALM.Workspace}, Ptr{Cfloat}, Ptr{Cfloat}),
            model.workspace,
            bmin != nothing ? pointer(bmin) : C_NULL,
            bmax != nothing ? pointer(bmax) : C_NULL,
        )
    end

    if length(kwargs) > 0
        settings = Settings(model; kwargs...)

        ccall(
            (:qpalm_update_settings, LIBQPALM_PATH),
            Cvoid,
            (Ptr{QPALM.Workspace}, Ptr{QPALM.Settings}),
            model.workspace, Ref(settings),
        )
    end
end

function copyto!(info::QPALM.Info, cinfo::QPALM.CInfo)
    info.iter = cinfo.iter
    info.iter_out = cinfo.iter_out
    info.status = QPALM.status_map[cinfo.status_val]
    info.status_val = cinfo.status_val
    info.pri_res_norm = cinfo.pri_res_norm
    info.dua_res_norm = cinfo.dua_res_norm
    info.dua2_res_norm = cinfo.dua2_res_norm
    info.objective = cinfo.objective
    info.dual_objective = cinfo.dual_objective
    info.setup_time = cinfo.setup_time
    info.solve_time = cinfo.solve_time
    info.run_time = cinfo.run_time

    return info
end

mutable struct Results
    x::Vector{Float64}
    y::Vector{Float64}
    info::QPALM.Info
    prim_inf_cert::Vector{Float64}
    dual_inf_cert::Vector{Float64}

    Results() = new(Float64[], Float64[], QPALM.Info(), Float64[], Float64[])
end

function warm_start!(model::QPALM.Model;
                     x_warm_start::Maybe{Vector{Float64}} = nothing,
                     y_warm_start::Maybe{Vector{Float64}} = nothing)
    if model.workspace == C_NULL
        error("QPALM Model not setup")
    end

    ccall(
        (:qpalm_warm_start, LIBQPALM_PATH),
        Cvoid,
        (Ptr{QPALM.Workspace}, Ptr{Float64}, Ptr{Float64}),
        model.workspace,
        x_warm_start != nothing ? pointer(x_warm_start) : C_NULL,
        y_warm_start != nothing ? pointer(y_warm_start) : C_NULL,
    )
end

function solve!(model::QPALM.Model, results::QPALM.Results=Results())::QPALM.Results
    ccall(
        (:qpalm_solve, LIBQPALM_PATH),
        Cvoid,
        (Ptr{QPALM.Workspace}, ),
        model.workspace
    )

    workspace = unsafe_load(model.workspace)

    info = unsafe_load(workspace.info)
    solution = unsafe_load(workspace.solution)
    data = unsafe_load(workspace.data)

    n = data.n
    m = data.m
    copyto!(results.info, info)

    resize!(results.x, n)
    resize!(results.y, m)
    resize!(results.prim_inf_cert, m)
    resize!(results.dual_inf_cert, n)

    has_solution = results.info.status in SOLUTION_PRESENT

    if has_solution
        # If solution exists, copy it
        unsafe_copyto!(pointer(results.x), solution.x, n)
        unsafe_copyto!(pointer(results.y), solution.y, m)
        fill!(results.prim_inf_cert, NaN)
        fill!(results.dual_inf_cert, NaN)
    else
        # else fill with NaN and return certificates of infeasibility
        fill!(results.x, NaN)
        fill!(results.y, NaN)
        if info.status == :Primal_infeasible || info.status == :Primal_infeasible_inaccurate
            unsafe_copyto!(pointer(results.prim_inf_cert), workspace.delta_y, m)
            fill!(results.dual_inf_cert, NaN)
        elseif info.status == :Dual_infeasible || info.status == :Dual_infeasible_inaccurate
            fill!(results.prim_inf_cert, NaN)
            unsafe_copyto!(pointer(results.dual_inf_cert), workspace.delta_x, n)
        else
            fill!(results.prim_inf_cert, NaN)
            fill!(results.dual_inf_cert, NaN)
        end
    end

    # An objective value is not calculated at this moment.
    # if results.info.status == :Non_convex
    #     results.info.obj_val = NaN
    # end

    return results
end

function cleanup!(model::QPALM.Model)
    ccall(
        (:qpalm_cleanup, LIBQPALM_PATH),
        Cvoid,
        (Ptr{QPALM.Workspace},),
        model.workspace
    )
end
