function refactor(d)
    N = size(d(1), 1)
    [namedtuple([k => d[k][i, :] for k = sort(keys(d))]) for i = 1:N]
end
_pmlfracs(a, N) = a[1:N, :]
_pmlfracs(a::Real, N) = fill(a, (N, 2))
_pmlfracs(a::AbstractVector, N) = stack([a[1:N], a[1:N]])


function setup(bbox, boundaries, sources, monitors, nres;
    approx_2D_mode=nothing,
    Ttrans=nothing, Tss=nothing, Tssmin=nothing,
    ϵ=1, μ=1, σ=0, m=0, γ=0, β=0,
    F=Float32,
    Courant=0.9,
    array=Array,
    pmlfracs=1,
    TEMP="",
)

    bbox, nres, = F.((bbox, nres,))
    rulers = makemesh(ϵ, bbox, nres)

    deltas = diff.(rulers)
    N = length(rulers)
    sz = Tuple(length.(rulers))

    if !isnothing(approx_2D_mode)
        approx_2D_mode = Symbol(approx_2D_mode)
    end
    pmlfracs = _pmlfracs(pmlfracs, N)

    geometry = OrderedDict()
    println("preprocessing geometry...")
    for (k, v) = pairs((; ϵ, μ, σ, m, γ, β))
        if isa(v, Number)
            v = F(v)
            if k ∈ (:σ, :m)
                geometry[k] = [fill(v, sz)]
            else
                geometry[k] = v
            end
        else
            geometry[k] = tensorinv(v, rulers)
            if k == :ϵ
                haspec = any(>=(PECVAL), last.(ϵ[1:end-1]))
                if !haspec
                    println("no PEC regions found in geometry")
                    geometry[:invϵ] = tensorinv(v, rulers; tensor=true, inv=true)
                else
                    geometry[:invϵ] = map(geometry[:ϵ][1]) do a
                        1 ./ a
                    end
                end
            elseif k == :μ
                geometry[:invμ] = map(geometry[:μ][1]) do a
                    1 ./ a
                end
            end
        end
    end

    a = geometry.ϵ[1]
    # a[a.>100]
    ϵmin, ϵmax = extrema(abs.(a))
    μmin, μmax = extrema(abs.(geometry.μ))
    nmax = sqrt(ϵmax * μmax)
    nmin = sqrt(ϵmin * μmin)

    dt = nmin / √(sum(v -> 1 / minimum(v)^2, deltas)) * Courant
    dt = 1 / ceil(1 / dt)

    maxdeltas = maximum.(deltas)

    v = 0.15 / dt
    δ = -3log(1e-4) / nmin / 2 / (2v) |> F
    σpml = ϵmin * v
    mpml = μmin * v

    pml_depths = max.(δ * pmlfracs, maxdeltas)
    pml_depths = ceil.(pml_depths ./ maxdeltas) .* maxdeltas
    @show σpml
    @show pml_depths

    geometry_names = (:ϵ, :μ, :σ, :m)
    if N == 1
        field_names = (:Ez, :Hy)
    elseif N == 2
        if approx_2D_mode == :TM
            Enames = (:Ez,)
            Hnames = (:Hx, :Hy)
            field_names = (:Ez, :Hx, :Hy)
        elseif approx_2D_mode == :TE
            Enames = (:Ex, :Ey)
            Hnames = (:Hz,)
            field_names = (:Ex, :Ey, :Hz)
        end
    else
        Enames = (:Ex, :Ey, :Ez)
        Hnames = (:Hx, :Hy, :Hz)
        field_names = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)
    end

    nodes = fill(:U, N, 2)
    db = Any[PML(j * i, pml_depths[i, int(j / 2 + 1.5)], σpml, mpml) for i = 1:N, j = (-1, 1)]

    nothingarray = () -> Array{Any,2}(fill(nothing, N, 2))
    zeroarray = () -> zeros(Int, N, 2)

    offsets = OrderedDict()
    sizes = OrderedDict([
        :default => collect(sz),
        (field_names .=> (collect(sz),))...
    ])
    boundvals = OrderedDict([
        :default => zeroarray(),
        (field_names .=> (nothingarray(),))...
    ])
    # 
    padvals = OrderedDict([
        :default => nothingarray(),
        (geometry_names .=> (nothingarray(),))...
    ])
    padamts = OrderedDict([
        :default => zeroarray(),
        (geometry_names .=> (zeroarray(),))...
    ])

    is_field_on_lb = OrderedDict([k => zeros(Int, N) for k = field_names])
    is_field_on_ub = OrderedDict([k => zeros(Int, N) for k = field_names])

    for b = boundaries
        for i = b.dims
            if typeof(b) == Periodic
                db[i, :] = [Periodic(-abs(i)), Periodic(abs(i))]
            else
                if i > 0
                    db[i, 2] = b
                else
                    db[abs(i), 1] = b

                end
            end
        end
    end

    for i = 1:N
        for j = 1:2
            b = db[i, j]
            t = typeof(b)
            if t == PML
            elseif t == PEC
                nodes[i, j] = :E
            elseif t == PMC
                nodes[i, j] = :H
            end
        end
    end

    for i = 1:N
        if nodes[i, :] == [:U, :E]

            nodes[i, 1] = :H
        elseif nodes[i, :] == [:U, :H]

            nodes[i, 1] = :E
        elseif nodes[i, :] == [:E, :U]

            nodes[i, 2] = :H
        elseif nodes[i, :] == [:H, :U]

            nodes[i, 2] = :E
        elseif nodes[i, :] == [:U, :U]

            nodes[i, :] = [:E, :H]
        elseif nodes[i, :] == [:U, :U]

            nodes[i, :] = [:E, :H]
        elseif nodes[i, :] == [:E, :E]

        elseif nodes[i, :] == [:H, :H]

        end

    end


    for i = 1:N
        select = i .== 1:N
        xyz = para = perp = [:x, :y, :z]
        perp = [popat!(para, i)]
        for j = 1:2
            b = db[i, j]
            if isa(b, PML)
                npml = round(b.d / maxdeltas[i])
                v = maxdeltas[i] * (1:npml)
                if j == 1
                    # pushfirst!(deltas[i], fill(maxdeltas[i], npml)...)
                    prepend!(rulers[i], rulers[i][1] - v)
                else
                    # push!(deltas[i], fill(maxdeltas[i], npml)...)
                    append!(rulers[i], rulers[i][end] + v)
                end

                padamts.default[i, j] = npml
                for k = geometry_names
                    padamts[k][i, j] = npml
                end

                for k = (:ϵ, :μ)
                    padvals[k][i, j] = :replicate
                end
                padvals[:σ][i, j] = TanhRamp(b.σ)
                padvals[:m][i, j] = TanhRamp(b.m)

                sizes.default[i] += npml
                for k = keys(sizes)
                    sizes[k][i] += npml
                end
            end

            f = nodes[i, j]
            for k = field_names
                q = startswith(String(k), String(f))
                if (q ? k[2] in para : k[2] in perp)
                    if isa(b, Periodic)
                        boundvals[k][i, j] = :periodic
                    else
                        boundvals[k][i, j] = 0
                    end
                end
            end
        end
    end


    for (k, v) = pairs(boundvals)
        is_field_on_lb[k] = !isnothing.(v[:, 1])
        is_field_on_ub[k] = !isnothing.(v[:, 2])
    end

    add_current_keys!(sizes)

    if N == 1
    elseif N == 3
        u0 = dict([k => zeros(F, Tuple(sizes[k])) for k = (:Ex, :Ey, :Ez, :Hx, :Hy, :Hz)])
    else
        if approx_2D_mode == :TM
            u0 = dict([k => zeros(F, Tuple(sizes[k])) for k = (:Ez, :Hx, :Hy)])
        else
            u0 = dict([k => zeros(F, Tuple(sizes[k])) for k = (:Ex, :Ey, :Hz)])
        end
    end
    add_current_keys!(u0)
    u0 = groupkeys(u0)

    offsets = OrderedDict([
        begin
            # xyz = f[2]
            # terminations = zip(is_field_on_lb[f], is_field_on_ub[f])
            # g = Symbol("$(k)$xyz$xyz")
            l = -is_field_on_lb[f] / 2
            f => [l l]
        end for f = field_names
    ])
    add_current_keys!(offsets)

    diffpadvals = kmap(a -> reverse(a, dims=2), boundvals)
    # diffpadvals = refactor(diffpadvals)

    rulers = OrderedDict([
        :default => rulers,
        begin
            d = map(rulers) do v
                vcat([2v[1] - v[2]], v, [2v[end] - v[end-1]])
            end
            map(pairs(offsets)) do (k, v)
                start = 2 + v[:, 1]
                stop = 1 + length.(rulers) + v[:, 2]
                k => getindexf.(d, range.(start, stop))
            end
        end...
    ])
    deltas = kmap(rulers) do vs
        diff.(vs)
    end
    dx = F(1 / nres / nmax)
    sizes = kmap(Tuple, sizes)
    grid = (; rulers, deltas, sizes, offsets, boundvals, diffpadvals, padvals, padamts, F, field_names, dx)

    println("making sources...")
    mode_solutions = []
    source_instances = SourceInstance.(sources, (grid,), (ϵ,), (TEMP,), (mode_solutions,),)
    println("making monitors...")
    monitor_instances = MonitorInstance.(monitors, (grid,), (ϵ,), (TEMP,), (mode_solutions,),)
    ϵeff = nothing

    if N == 2
        # geometry[:ϵ] = downsample(_geometry.ϵ, int(deltas / dl))
    end

    if !isa(Ttrans, Real)
        if endswith(string(Ttrans), "x")
            Ttrans = parse(F, Ttrans[1:end-1])
        else
            Ttrans = 1
        end
        Ttrans *= sum(L * nmax)
    end

    if Tss == nothing
        if isnothing(Tssmin)
            Tssmin = 120nmax / nres / N
        end
        v = reduce(vcat, wavelengths.(monitor_instances))
        v = Base.round.(v, digits=3)
        v = v |> Set |> collect |> sort |> reverse
        @show v
        if length(v) == 1
            Tss = 1
        else
            Tss = 1 / minimum(diff([0, (1 ./ v)...]))
            # Tss = ceil(100 / T) * T
        end
    end
    if !isnothing(Tssmin)
        @show Tssmin
        Tss *= Base.ceil(Int, Tssmin / Tss)
    end

    @show Ttrans, Tss
    Ttrans, Tss = convert.(F, (Ttrans, Tss))
    global prob = (;
                      grid,
                      source_instances, monitor_instances, field_names, approx_2D_mode, Courant,
                      Ttrans, Tss,
                      geometry, nmax, nmin, ϵeff,
                      u0, dt, array,) |> pairs |> OrderedDict

    _gpu = x -> gpu(array, x)
    if array == Array
        println("using CPU backend.")
    else
        println("using GPU backend.")
    end
    for k in (:u0, :_geometry, :geometry, :source_instances, :monitor_instances)
        prob[k] = _gpu(prob[k],)
    end
    for k = [:diffdeltas]
        prob.grid[k] = _gpu.(prob.grid[k])
    end

    geometry = pad_geometry(geometry, padvals, padamts)

    p = fmap(array, p, AbstractArray{<:Number})

    prob
end
update = update
setup = setup