function refactor(d)
    N = size(d(1), 1)
    [namedtuple([k => d[k][i, :] for k = sort(keys(d))]) for i = 1:N]
end
_pmlfracs(a, N) = a[1:N, :]
_pmlfracs(a::Real, N) = fill(a, (N, 2))
_pmlfracs(a::AbstractVector, N) = stack([a[1:N], a[1:N]])


function setup(λ, bbox, nres, boundaries, sources, monitors, canvases=[];
    approx_2D_mode=nothing, z=nothing,
    Ttrans=nothing, Tss=nothing, Tssmin=nothing,
    ϵ=1, μ=1, σ=0, m=0, γ=0, β=0,
    F=Float32,
    relcourant=0.9,
    array=Array,
    pmlfracs=1,
    TEMP="",
)
    println("setting up simulation...")

    bbox /= λ
    if !isnothing(z)
        z /= λ
    end

    bbox, nres, relcourant = F.((bbox, nres, relcourant))

    T = Float32
    ϵ = map(ϵ) do (m, v)
        isnothing(m) ? m : m |> Scale(1 / λ, 1 / λ, 1 / λ), v
    end
    rulers = makemesh([(m, sqrt(v) |> T) for (m, v) = ϵ], bbox |> T, nres |> T) |> F

    deltas = diff.(rulers)
    @debug extrema.(deltas)
    @debug deltas[end]
    N = length(rulers)
    sz = Tuple(length.(rulers) - 1)
    # error("stop here")
    if !isnothing(approx_2D_mode)
        approx_2D_mode = Symbol(approx_2D_mode)
    end
    pmlfracs = _pmlfracs(pmlfracs, N)

    geometry = OrderedDict()
    println("meshing geometry - can take few minutes...")
    for (k, v) = pairs((; ϵ, μ, σ, m, γ, β))
        if isa(v, Number)
            v = F(v)
            if k ∈ (:σ, :m)
                geometry[k] = [fill(v, sz)]
            else
                geometry[k] = v
                if k == :μ
                    geometry[:invμ] = 1 / v
                end
            end
        else
            if k == :ϵ
                haspec = any(>=(PECVAL), last.(ϵ[1:end-1]))
                if !haspec
                    # println("no PEC regions found in geometry")
                    @time geometry[:invϵ] = supersamplemesh(v, rulers; tensor=true, inv=true, z, big=!AUTODIFF())
                    # geometry[:invϵ] = [map(geometry[:ϵ][1]) do a
                    #     1 ./ a
                    # end]
                    geometry[:ϵ] = [1 ./ mean(diag(geometry[:invϵ]))]
                else
                    geometry[:ϵ] = supersamplemesh(v, rulers; z)
                    geometry[:invϵ] = [1 ./ geometry[:ϵ][1]]
                end
            else
                geometry[k] = supersamplemesh(v, rulers; z)
                if k == :μ
                    geometry[:invμ] = [map(geometry[:μ][1]) do a
                        1 ./ a
                    end]
                end
            end
        end
    end
    geometry_names = (:ϵ, :μ, :σ, :m, :invϵ, :invμ)
    for k = geometry_names
        @debug k extrema(geometry[k][1])
    end

    a = geometry.ϵ[1]
    # a[a.>100]
    ϵmin, ϵmax = extrema(abs.(a))
    μmin, μmax = extrema(abs.(geometry.μ))
    nmax = sqrt(ϵmax * μmax)
    nmin = sqrt(ϵmin * μmin)

    dt = nmin / √(sum(v -> 1 / minimum(v)^2, deltas)) * relcourant
    dt = 1 / ceil(1 / dt) |> F


    v = 0.15 / dt |> F
    δ = -2log(1.0f-4) / nmin / 2 / (2v) |> F
    σpml = ϵmin * v
    mpml = μmin * v

    pml_depths = δ * pmlfracs
    @debug σpml pml_depths

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
        map(field_names) do k
            k => collect(sz)
        end...
    ])
    boundvals = OrderedDict([
        # :default => zeroarray(),
        map(field_names) do k
            k => nothingarray()
        end...
    ])
    # 
    padvals = OrderedDict([
        :default => nothingarray(),
        map(geometry_names) do k
            k => nothingarray()
        end...
    ])
    padamts = OrderedDict([
        :default => zeroarray(),
        map(geometry_names) do k
            k => zeroarray()
        end...
    ])

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
                d = j == 1 ? deltas[i][1] : deltas[i][end]
                npml = round(b.d / d)
                npml = max(npml, 1)
                v = d * (1:npml)
                if j == 1
                    # pushfirst!(deltas[i], fill(maxdeltas[i], npml)...)
                    prepend!(rulers[i], rulers[i][1] - reverse(v))
                else
                    # push!(deltas[i], fill(maxdeltas[i], npml)...)
                    append!(rulers[i], rulers[i][end] + v)
                end

                padamts.default[i, j] = npml
                for k = geometry_names
                    padamts[k][i, j] = npml
                end

                for k = (:ϵ, :μ, :invϵ, :invμ)
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

    edges = OrderedDict()
    for (k, v) = pairs(boundvals)
        edges[k] = !isnothing.(v)
    end
    edges[:default] = zeros(Bool, N, 2)

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
    all_field_names = keys(u0)
    u0 = groupkeys(u0)

    add_current_keys!(edges)
    offsets = OrderedDict([
        begin
            # xyz = f[2]
            # terminations = zip(is_field_on_lb[f], is_field_on_ub[f])
            # g = Symbol("$(k)$xyz$xyz")
            if f != :default
                l = -edges[f][:, 1] / 2 + 0.25
                v = [l l] |> F
            end
            f => v
        end for (f, v) = pairs(edges)])

    diffpadvals = vmap(a -> reverse(a, dims=2), boundvals)

    rulers = OrderedDict([
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
    deltas = vmap(rulers) do vs
        diff.(vs)
    end
    diffdeltas = OrderedDict([
        k => begin
            vs = centroids.(vs)
            N = length(vs)
            map(1:N, vs, eachrow(edges[k])) do i, v, (l, r)
                s = i .== (1:N)
                v = diff(v)
                if !l
                    pushfirst!(v, first(v))
                end
                if !r
                    push!(v, last(v))
                end
                reshape(v, (.!s + length(v) * s)...)
            end
        end for (k, vs) = pairs(rulers)])

    dx = F(1 / nres / nmax)
    sizes = vmap(Tuple, sizes)
    grid = (; rulers, deltas, sizes, edges, offsets, boundvals, diffpadvals, diffdeltas, padvals, padamts, F, N, field_names, all_field_names, dx) |> pairs |> OrderedDict
    geometry = pad_geometry(geometry, padvals, padamts) |> pairs |> OrderedDict


    println("making sources...")
    mode_solutions = []
    source_instances = SourceInstance.(sources, λ, (grid,), (ϵ,), (TEMP,); z, mode_solutions,)
    println("making monitors...")
    monitor_instances = MonitorInstance.(monitors, λ, (grid,), (ϵ,), (TEMP,); z, mode_solutions,)
    canvas_instances = CanvasInstance.(canvases, λ, (grid,); z)

    if N == 2
        # geometry[:ϵ] = downsample(_geometry.ϵ, int(deltas / dl))
    end
    dimensions = last.(rulers.default) - first.(rulers.default)
    if !isa(Ttrans, Real)
        if endswith(string(Ttrans), "x")
            Ttrans = parse(F, Ttrans[1:end-1])
        else
            Ttrans = 1
        end
        Ttrans *= sum(dimensions * nmax) + 4
    end

    if Tss == nothing
        if isnothing(Tssmin)
            Tssmin = max(10, 50nmax / nres / N)
        end
        v = reduce(vcat, wavelengths.(monitor_instances))
        v = Base.round.(v, digits=3)
        v = v |> Set |> collect |> sort |> reverse
        if length(v) == 1
            Tss = 1
        else
            Tss = 1 / minimum(diff([0, (1 ./ v)...]))
            # Tss = ceil(100 / T) * T
        end
    end
    if !isnothing(Tssmin)
        Tss *= Base.ceil(Int, Tssmin / Tss)
    end
    Ttrans, Tss = convert.(F, (Ttrans, Tss))
    Tss = dt * round(Tss / dt)
    Ttrans = dt * round(Ttrans / dt)

    T = Ttrans + Tss
    nt = round(Int, T / dt)
    nv = prod(sizes.default)
    load = nt * nv

    global prob = (;
                      grid, λ,
                      source_instances, monitor_instances, canvas_instances,
                      field_names, approx_2D_mode,
                      Ttrans, Tss,
                      geometry, nmax, nmin,
                      u0, dt, array,) |> pairs |> OrderedDict
    if array == Array
        backend = :CPU
    else
        backend = :GPU
        for k in (:u0, :geometry, :source_instances, :monitor_instances)
            prob[k] = gpu(array, prob[k],)
        end
        for k = [:diffdeltas]
            prob.grid[k] = gpu(array, prob.grid[k])
        end
    end
    println(BREAK)
    println()

    println("simulation setup complete")
    println()
    # println("Courant number: $Courant")
    println("backend: $backend")
    println("float: $F")
    println()

    println("original size: $sz")
    println("padded size: $(sizes.default)")
    println("cell count: $(nv|>disp)")
    println()

    println("transient time: $(Ttrans|>disp) periods")
    println("steady state time: $(Tss|>disp) periods")
    println("total time: $(T|>disp) periods")
    println("time steps: $(nt|>disp)")

    println()
    println("computation load: $(load|>disp) cell-steps")
    println(DBREAK)
    println("")
    prob
end
update = update
setup = setup