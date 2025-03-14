function bell(t, dt, u=nothing)
    ignore_derivatives() do

        # fmap(u) do x
        #     isnan(x) && error("NaN detected in field")
        # end

        if floor(t + 0.001) > floor(t - dt + 0.001)
            !AUTODIFF() && t > 0.1 && println("period $(int(t)),  $(timepassed()|>disp) s")
        end
        if !haskey(ENV, "t0")
            ENV["t0"] = time()
        end
    end
end

function f1(((u,), p, (dt, diffdeltas, diffpadvals, source_instances,)), t)
    bell(t, dt, u)
    # @time u = update(u, p, t, dt, diffdeltas, diffpadvals, source_instances)
    u = update(u, p, t, dt, diffdeltas, diffpadvals, source_instances;)
    ((u,), p, (dt, diffdeltas, diffpadvals, source_instances,))
end

function f2(((u, um), p, (dt, diffdeltas, diffpadvals, source_instances,), (t0, T, monitor_instances)), t)
    bell(t, dt, u)
    # @time u = update(u, p, t, dt, diffdeltas, diffpadvals, source_instances;)
    u = update(u, p, t, dt, diffdeltas, diffpadvals, source_instances;)
    ks = @ignore_derivatives [keys(u.E)..., keys(u.H)...]
    um += [
        begin
            um = namedtuple(ks .=> field.((u,), ks, (m,)))
            [
                begin
                    c = dt / T * cispi(-2(t - t0) / λ)
                    @nograd c
                    c * um
                end for λ = wavelengths(m)
            ]
        end for m = monitor_instances]
    ((u, um), p, (dt, diffdeltas, diffpadvals, source_instances,), (t0, T, monitor_instances))
end

function solve(prob, models=nothing;
    save_memory=false, ulims=(-3, 3), framerate=nothing, showfield=:Hz, path="", verbose=true,
    kwargs...)
    @unpack approx_2D_mode, dt, u0, geometry, source_instances, monitor_instances, Ttrans, canvas_instances, Tss, grid, array = prob

    @nograd dt, u0, source_instances, monitor_instances, Ttrans, Tss, grid
    @unpack diffdeltas, diffpadvals, F, N, = grid

    for c = canvas_instances
        apply_canvas!(c, geometry, models)
    end
    p = geometry

    TEMP = joinpath(path, "temp")
    FIELDS = joinpath(path, "temp", "fields")
    mkpath(FIELDS)
    mkpath(TEMP)

    if !isnothing(framerate)
        npzwrite(joinpath(TEMP, "g.npy"), prob._geometry.ϵ)
    end


    durations = [Ttrans, Tss]
    T = cumsum(durations)
    us0 = (u0,)
    init = (us0, p, (dt, diffdeltas, diffpadvals, source_instances,))

    ts = range(0, T[1] - dt, int(durations[1] / dt))
    @nograd ts

    @ignore_derivatives delete!(ENV, "t0")

    verbose && println("propagating transient fields...")
    timepassed()
    if save_memory
        (u,), = adjoint_reduce(f1, ts, init, ulims)
    else
        (u,), = reduce(ts; init) do us, t

            ignore() do
                if !isnothing(framerate) && t > 0
                    if t % (1 / framerate) < dt
                        (u,), = us
                        a = u(showfield)
                        npzwrite(joinpath(FIELDS, "$t.npy"), a)

                        # CairoMakie.save(joinpath(FIELDS, "$t.png"), quickie(a, g; monitor_instances, source_instances, ulims),)
                        # quickie(a, g; monitor_instances, source_instances)
                    end
                end
            end
            f1(us, t;)
        end
    end

    if array == Array
        fmap(u) do x
            isnan(x) && error("NaN detected in field")
        end
    end

    ts = range(T[1], T[2] - dt, int(durations[2] / dt))
    init = ((u, 0), p, (dt, diffdeltas, diffpadvals, source_instances,), (T[1], durations[2], monitor_instances,))
    @nograd ts


    verbose && println("accumulating dft fields...")
    if save_memory
        (u, um), = adjoint_reduce(f2, ts, init, ulims)
    else
        (u, um), = reduce(f2, ts; init)
    end
    # return (u.E.Ex + u.H.Hz + u.E.Ey) .|> abs |> sum

    ignore_derivatives() do
        t0 = parse(Float64, ENV["t0"])
        println("simulation done in $(time() - t0|>disp) s (includes JIT compilation).")
    end

    # return sum(abs, um[1][1].Ex + um[1][1].Ey + um[1][1].Hz)
    ulims = 0
    # @assert all([all(!isnan, a) for a = u])

    # volume(cpu(prob.monitor_instances[2].(1)) + 0.001ϵ(1) |> cpu) |> display
    # extrema(cpu(prob.monitor_instances[2].λmodes(1)(1).Ey))
    # extrema(cpu(prob.monitor_instances[2].λmodes(1)(1).Hx))
    # volume(cpu(abs.(prob.source_instances[1].sigmodes(1)[2].Jy))) |> display
    # for i = 1:2
    #     for k = (:Ey, :Hx)
    #         prob.monitor_instances[i]._λmodes(1)[1](k) |> extrema |> println
    #     end
    # end
    # global a = um
    # conj(a[1][1].Ey) .* a[1][1].Hx |> sum |> println
    # conj(a[2][1].Ey) .* a[2][1].Hx |> sum |> println
    # error()

    @nograd monitor_instances
    v = map(um, monitor_instances) do um, m
        @nograd m
        map(um, wavelengths(m)) do um, λ
            um, m = cpu(um), cpu(m)
            um = localframe(um, m; approx_2D_mode)
            ap = am = nothing
            if !isnothing(m.λmodes)
                modes = m.λmodes[λ]
                _modes = m._λmodes[λ]
                ap = inner.(modes, (um,), (m.plane_deltas,))
                am = inner.(_modes, (um,), (m.plane_deltas,))
            end
            P = real(inner(um, um, m.plane_deltas))
            um, ap, am, P
        end
    end

    # extrema(prob.monitor_instances[1].λmodes(1)(1).Hx)
    um = [[v[1] for v = v] for v in v]
    ap = [[v[2] for v = v] for v in v]
    am = [[v[3] for v = v] for v in v]
    P = [[v[4] for v = v] for v in v]

    # nm = length(monitor_instances)
    # nλ = length(wavelengths(monitor_instances[1]))
    # um = [[v[j][i][1] for j = 1:nλ] for i = 1:nm]
    # ap = [[v[j][i][2] for j = 1:nλ] for i = 1:nm]

    # am = [[v[j][i][3] for j = 1:nλ] for i = 1:nm]

    return Solution(u, p, um, P, ap, am)
end


struct Solution
    u
    p
    um
    P
    ap
    am
end
@functor Solution

function (s::Solution)(k, m, w=1, mn=0)
    @unpack u, P, um, ap, am = s
    if k == "a+"
        return s.ap[m][w][mn+1]
    elseif k == "a-"
        return s.am[m][w][mn+1]
    elseif k == "P_TE"
        return flux(um[m][w], :TE)
    elseif k == "P_TM"
        return flux(um[m][w], :TM)
    elseif k == "P"
        return P[m][w]
    elseif k == "um"
        return um[m][w]
    end
end
# heatmap(___p.invϵ[1, 1])

function getsparams(s, iλ=1, imode=1)
    @unpack u, ulims, um, ap, am = s
    # nλ = length(ap[1])
    np = length(ap)
    # ifelse.(1:np .== 1, am[j], ap[j]) 
    map(1:np) do i
        ap[i][iλ][imode]
    end / am[1][iλ][imode]
end