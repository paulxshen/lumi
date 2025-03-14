function picrun(path)
    PROB = joinpath(path, "problem.json")
    prob = JSON.parse(read(open(PROB), String); dicttype=OrderedDict)
    gpu_backend = prob["gpu_backend"]
    array = if isnothing(gpu_backend)
        Array
    else
        cu
    end

    # function picrun(path, array=Array; kw...)
    Random.seed!(1)
    PROB = joinpath(path, "problem.json")
    SOL = joinpath(path, "solution.json")
    TEMP = joinpath(path, "temp")

    io = open(PROB)
    s = read(io, String)
    prob = JSON.parse(s; dicttype=OrderedDict)
    # merge!(prob, kw)
    # for (k, v) = pairs(kw)
    #     prob[string(k)] = v
    # end
    @unpack name, N, approx_2D_mode, dtype, center_wavelength, runs, ports, study, zmin, zmax, zcenter, magic, framerate, layer_stack, materials, Ttrans, Tss, bbox, epdefault, nres = prob
    if nres > 8
        error("nres must be <= 8 in free version")
    end

    if study == "inverse_design"
        @unpack canvases, lsolid, lvoid, targets, weights, eta, iters, restart, save_memory, design_config, stoploss = prob
    end

    F = Float32
    dtype = lowercase(dtype)
    if contains(dtype, "16") && contains(dtype, "bf")
        F = BFloat16
        println("BFloat16 selected. make sure your GPU supports it. otherwise will be emulated and very slow.")

    elseif contains(dtype, "16")
        F = Float16
        println("Float16 selected. make sure your cpu or GPU supports it. otherwise will be emulated and very slow.")
    end

    bbox = stack(bbox)
    @show λ = F(center_wavelength)

    layer_stack = sort(collect(pairs(layer_stack)), by=kv -> kv[2].mesh_order) |> OrderedDict
    ks = keys(layer_stack)
    fns = readdir(joinpath(path, "bodies"), join=true)
    sort!(fns)
    @debug fns
    meshes = getfield.(GeoIO.load.(fns, numbertype=Float32), :domain)
    eps = [materials(string(split(basename(fn), "_")[2])).epsilon |> F for fn = fns]
    meps = zip(meshes, eps)
    epdefault = F(epdefault)

    models = nothing
    z = nothing
    if N == 2
        # midplane = Plane((0, 0, zcenter), (0, 0, 1))
        # meps = map(meps) do (m, ϵ)
        #     m ∩ midplane, ϵ
        # end
        z = mean(bbox[3, :])
        bbox3 = bbox
        bbox = bbox[1:2, :]
    end
    meps = collect(Tuple{Any,Any}, meps)
    push!(meps, (nothing, epdefault))
    ϵ = meps

    if study == "inverse_design"
        ENV["AUTODIFF"] = 1
        targets = fmap(F, targets)
        # targets = sortkeys(targets)
        if isfile(SOL)
            sol = open(SOL, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
    end

    boundaries = [] # unspecified boundaries default to PML

    runs = [SortedDict([k => isa(v, AbstractDict) ? SortedDict(v) : v for (k, v) = pairs(run)]) for run in runs]
    runs_sources = [
        begin
            sources = []
            for (port, sig) = SortedDict(run.sources) |> pairs
                @unpack center, frame, dimensions, wavelength_mode_numbers = sig
                center = center[1:N]
                dimensions = dimensions[1:N-1]
                frame = stack(frame)
                if N == 2
                    frame = [frame[1:2, 1] frame[1:2, 3]]
                end
                λmodenums = SortedDict([(F(λ)) => v for (λ, v) in pairs(wavelength_mode_numbers)])

                push!(sources, Source(center, dimensions, frame; λmodenums, port, wavelength_mode_numbers, label="s$(string(port)[2:end])"))
            end
            sources
        end for run in runs
    ]
    # sort!(runs_sources, by=x -> x.label)

    runs_monitors = [[
        begin
            @unpack center, frame, dimensions, wavelength_mode_numbers = m
            center = center[1:N]
            dimensions = dimensions[1:N-1]
            frame = stack(frame)
            if N == 2
                frame = [frame[1:2, 1] frame[1:2, 3]]
            end
            λmodenums = SortedDict([F(λ) => v for (λ, v) in pairs(m.wavelength_mode_numbers)])

            Monitor(center, dimensions, frame; λmodenums, port, label="m$(string(port)[2:end])", wavelength_mode_numbers)
        end for (port, m) = SortedDict(run.monitors) |> pairs] for run in runs]
    if study == "inverse_design"
        canvases = map(enumerate(canvases)) do (i, c)
            @unpack symmetries, swaps = c
            if !isnothing(sol) && !restart
                println("loading saved design...")
                params = sol.params[i]
            else
                params = nothing
            end
            swaps = kvmap(swaps) do k, v
                Symbol(togreek(k)) => v
            end
            Canvas(swaps, stack(c.bbox), lvoid, lsolid, symmetries; params)
        end
    else
        canvases = []
    end
    probs =
        [
            begin
                setup(λ, bbox, nres, boundaries, sources, monitors, canvases;
                    pmlfracs=[1, 1, 0.2], approx_2D_mode, z, array,
                    F, ϵ, TEMP, Ttrans, Tss)
            end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
        ]

    # error("not implemented")
    t0 = time()
    println("compiling simulation code...")
    if study == "sparams"
        res = calc_sparams(probs; verbose=true)
        @unpack S, T, sols = res
        plotsim(probs[1] |> cpu, sols[1] |> cpu, ; path=joinpath(path, "sim.png"))

        sol = (;
            sparam_family(S)...,
            total_power_T=T,
            total_power_dB=todB(T),
            path, study)
        println("modal T-params: ")
        println(json(sol.T, 4))
        println("total power T-params: ")
        println(json(T, 4))
        open(SOL, "w") do f
            write(f, json(cpu(sol)))
        end
    elseif study == "inverse_design"
        PLOT = joinpath(path, "sim.png")
        if N == 3
            if magic != "summersale"
                error("3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")
            end
        end
        prob = probs[1]
        model = prob.canvas_instances[1].model
        # minchange = 0.001
        opt = AreaChangeOptimiser(
            model; opt=Optimisers.Momentum(1, 0.7),
        )
        opt_state = Flux.setup(opt, model)
        # error("not implemented")
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        println("")
        # iters = 1
        for i = 1:iters
            println("[$i]  ")
            stop = i == iters
            function f(model)
                models = [model]
                r = calc_sparams(probs, models)
                global S = r.S
                global T = r.T
                global sols = r.sols
                @debug targets
                l = 0
                for k = keys(targets)
                    y = targets[k]
                    err = (x, y) -> abs(x - y)
                    if :phasediff == k
                        yhat = namedtuple([
                            _λ => namedtuple([
                                ps =>
                                    begin
                                        ks = ignore_derivatives() do
                                            ps = split(string(ps), ",")
                                            ks = keys(S[_λ])
                                            [ks[findfirst(k -> startswith(string(k), p), ks)] for p = ps]
                                        end
                                        s1, s2 = [S[_λ][k] for k in ks]
                                        angle(s1 / s2)
                                    end
                                for ps = keys(targets[k][_λ])])
                            for _λ = keys(targets[k])])
                        err = (x, y) -> angle(cis(x) / cis(y))
                        yhat = flatten(yhat)
                        y = flatten(y)
                        Z = length(y) * convert.(F, 2π)
                    else
                        yhat = if "sparams" == string(k)
                            S
                        elseif "tparams" == string(k)
                            fmap(abs2, S)
                        end

                        yhat = [[yhat(_λ)(k) for k = keys(y[_λ])] for _λ = keys(y)]
                        yhat = flatten(yhat)
                        y = flatten(y)
                        Z = sum(abs, y)
                    end
                    v = err.(yhat, y) / Z
                    println("$(k) losses: $v ")
                    l += sum(v)
                end
                println("modified total loss: $l\n")
                println()
                l
            end

            if iters == 1
                f(model)
                break
            else
                @time l, (dldm,) = Flux.withgradient(f, model)
            end
            @assert !isnothing(dldm)
            if !isnothing(stoploss) && l < stoploss
                println("\nLoss below threshold, stopping optimization.")
                stop = true
            end
            if true# i == 1 || i % 2 == 0 || stop
                println("saving checkpoint...")
                CKPT = joinpath(path, "checkpoints", replace(string(now()), ':' => '_', '.' => '_'))

                mkpath(CKPT)
                models = [model]
                global S, T, sols
                sol = (;
                    sparam_family(S)...,
                    total_power_T=T,
                    total_power_dB=todB(T),
                    params=map(prob.canvas_instances) do c
                        c.model.p |> cpu
                    end)
                open(SOL, "w") do f
                    write(f, json(cpu(sol)))
                end

                open(joinpath(CKPT, "solution.json"), "w") do f
                    write(f, json(cpu(sol)))
                end
                println("T-params: ")
                println(json(sol.T, 4))
                plotsim(probs[1] |> cpu, sols[1] |> cpu, ; path=PLOT)

            end
            if stop
                break
            end
            # opt.minchange = i == 1 ? 0.1 : 0.001
            opt.minchange = max(0.001, 0.02l^2)
            # opt.maxchange = maxchange * (1 + 2l)
            Jello.update_loss!(opt, l)
            Flux.update!(opt_state, model, dldm)# |> gpu)
            GC.gc(true)
            println(BREAK)
            println()
        end
        println("iteration done in $((time() - t0)|>disp) s")
    end
end