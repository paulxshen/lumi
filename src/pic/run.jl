# function picrun(path)
#     PROB = joinpath(path, "problem.json")
#     prob = JSON.parse(read(open(PROB), String); dicttype=OrderedDict)
#     gpu_backend = prob["gpu_backend"]
#     array = if isnothing(gpu_backend)
#         println("using CPU")
#         Array
#     else
#         println("using $gpu_backend")
#         cu
#     end
#     picrun(path, array)
# end

function picrun(path, array=Array; kw...)
    Random.seed!(1)
    ENV["autodiff"] = "0"
    println("setting up simulation...")
    PROB = joinpath(path, "problem.json")
    SOL = joinpath(path, "solution.json")
    TEMP = joinpath(path, "temp")

    io = open(PROB)
    s = read(io, String)
    prob = JSON.parse(s; dicttype=OrderedDict)
    # merge!(prob, kw)
    for (k, v) = pairs(kw)
        prob[string(k)] = v
    end
    @unpack name, N, approx_2D_mode, dtype, center_wavelength, runs, ports, study, zmin, zmax, zcenter, magic, framerate, layer_stack, materials, Ttrans, Tss, bbox, epdefault, nres = prob
    if study == "inverse_design"
        @unpack lsolid, lvoid, designs, targets, weights, eta, iters, restart, save_memory, design_config, stoploss = prob
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
    println("using $F")

    bbox = stack(bbox)

    @show λ = F(center_wavelength)

    layer_stack = sort(collect(pairs(layer_stack)), by=kv -> kv[2].mesh_order) |> OrderedDict
    ks = keys(layer_stack)
    fns = readdir(joinpath(path, "surfaces"), join=true)
    sort!(fns)
    @show fns
    global meshes = getfield.(GeoIO.load.(fns, numbertype=F), :domain) .|> (Scale(1 / λ, 1 / λ, 1 / λ,),)
    global eps = [materials(string(split(basename(fn), "_")[2])).epsilon |> F for fn = fns]
    meps = zip(meshes, eps)
    epdefault = F(epdefault)

    global models = nothing
    if N == 2
        midplane = Plane((0, 0, zcenter), (0, 0, 1))
        meps = map(meps) do (m, ϵ)
            m ∩ midplane, ϵ
        end
    end
    meps = collect(Tuple{Any,Any}, meps)
    push!(meps, (nothing, epdefault))
    ϵ = meps

    if study == "inverse_design"
        targets = fmap(F, targets)
        # targets = sortkeys(targets)
        if isfile(SOL)
            sol = open(SOL, "r") do f
                JSON.parse(f)
            end
        else
            sol = nothing
        end
        models = [
            begin
                @unpack bbox = design
                dimensions = bbox[2] - bbox[1]
                szd = Tuple(round.(Int, dimensions / dl)) # design region size
                symmetries = map(design.symmetries) do s
                    try
                        Int(s) + 1
                    catch
                        try
                            parse(Int, s) + 1
                        catch
                            s
                        end
                    end
                end

                frame = ϵ2
                # display(heatmap(frame))
                # error("not implemented")
                frame = frame .>= 0.99maximum(frame)
                # frame = nothing
                start = round((bbox[1] - lb) / dl + 1)
                b = Blob(szd; solid_frac=1, lsolid=lsolid / dl, lvoid=lvoid / dl, symmetries, F, frame, start, morph=false)
                display(heatmap(b.frame))

                if !isnothing(sol) && !restart
                    println("loading saved design...")
                    b.p .= sol.params[i]
                end
                b
            end for (i, design) = enumerate(designs)
        ]
    end

    boundaries = [] # unspecified boundaries default to PML

    λ = F(λ)
    global runs = [SortedDict([k => isa(v, AbstractDict) ? SortedDict(v) : v for (k, v) = pairs(run)]) for run in runs]
    global runs_sources = [
        begin
            sources = []
            for (port, sig) = SortedDict(run.sources) |> pairs
                @unpack center, frame, dimensions, wavelength_mode_numbers = sig
                center = center[1:N]
                dimensions = dimensions[1:N-1]
                frame = stack(frame)

                λmodenums = SortedDict([(F(λ)) => v for (λ, v) in pairs(wavelength_mode_numbers)])

                push!(sources, Source(center, dimensions, frame, λ; λmodenums, label="s$(string(port)[2:end])"))
            end
            sources
        end for run in runs
    ]
    # sort!(runs_sources, by=x -> x.label)

    global runs_monitors = [[
        begin
            @unpack center, frame, dimensions, wavelength_mode_numbers = m
            center = center[1:N]
            dimensions = dimensions[1:N-1]
            frame = stack(frame)

            λmodenums = SortedDict([F(λ) => v for (λ, v) in pairs(m.wavelength_mode_numbers)])

            Monitor(center, dimensions, frame, λ; λmodenums, label=port)
        end for (port, m) = SortedDict(run.monitors) |> pairs] for run in runs]

    global run_probs =
        [
            begin
                setup(bbox / λ, boundaries, sources, monitors, nres;
                    pmlfracs=[1, 1, 0.2], approx_2D_mode, array,
                    F, ϵ, TEMP, Ttrans, Tss)
            end for (i, (run, sources, monitors)) in enumerate(zip(runs, runs_sources, runs_monitors))
        ]

    # error("not implemented")
    t0 = time()
    println("compiling simulation code...")
    if study == "sparams"
        @unpack S, sols = calc_sparams(runs, run_probs, lb, dl;
            F, verbose=true, framerate, path)
        plotsols(sols, run_probs, path)
        sol = (; sparam_family(S)...,
            path, study)
        open(SOL, "w") do f
            write(f, json(cpu(sol)))
        end
    elseif study == "inverse_design"
        ENV["autodiff"] = "1"
        if length(lb) == 3
            if magic != "summersale"
                error("3D inverse design feature must be requested from Luminescent AI info@luminescentai.com")
            end
        end
        model = models[1]
        # minchange = 0.001
        # maxchange = max(4minchange, holesize(model) / prod(size(model)))
        global opt = AreaChangeOptimiser(
            model; opt=Optimisers.Momentum(1, 0.7),
        )
        opt_state = Flux.setup(opt, model)
        # error("not implemented")
        println("starting optimization... first iter will be slow due to adjoint compilation.")
        img = nothing
        println("")
        for i = 1:iters
            println("====($i)  ")
            stop = i == iters
            if :phase_shifter == first(keys(targets))
                @time l, (dldm,) = Flux.withgradient(model) do model
                    models = [model]
                    @unpack S, sols = make_pic_sim_problem(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg)#(1)(1)
                    k = keys(S) |> first
                    s = S[k][Symbol("o2@0,o1@0")]

                    @unpack S = make_pic_sim_problem(runs, run_probs, lb, dl,
                        designs, design_config, models;
                        F, img, alg, save_memory, perturb=:ϵ,)#(1)(1)
                    s_ = S[k][Symbol("o2@0,o1@0")]

                    T = abs2(s)
                    dϕ = angle(s_ / s)
                    println("T: $T, dϕ: $dϕ")
                    (exp(T - 1) * dϕ / π)
                    # T * dϕ / π
                end
            else
                function f(model)
                    models = [model]
                    res = calc_sparams(runs, run_probs, lb, dl,
                        designs, design_config, models, ;
                        F, img, alg, save_memory, materials, path, framerate)
                    # return res
                    @unpack S, sols = res
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

                            # global a1 = S, y
                            yhat = [[yhat(_λ)(k) for k = keys(y[_λ])] for _λ = keys(y)]
                            yhat = flatten(yhat)
                            # println("yhat: $yhat")
                            y = flatten(y)
                            # println("y: $y")
                            Z = sum(abs, y)
                        end
                        # _l = sum(,) do x
                        #     (x) + x^2
                        # end
                        v = err.(yhat, y) / Z
                        println("$(k) losses: $v ")
                        # v = v .* @ignore_derivatives softmax(20v) * length(v)
                        l += sum(v)
                    end
                    println("modified total loss: $l\n")
                    println()
                    l
                end

                # f(model)
                if iters == 1
                    f(model)
                    break
                else
                    @time global l, (dldm,) = Flux.withgradient(f, model)
                end
            end
            @assert !isnothing(dldm)
            if !isnothing(stoploss) && l < stoploss
                println("\nLoss below threshold, stopping optimization.")
                stop = true
            end
            if true# i == 1 || i % 2 == 0 || stop
                println("saving checkpoint...")
                ckptpath = joinpath(path, "checkpoints", replace(string(now()), ':' => '_', '.' => '_'))

                mkpath(ckptpath)
                models = [model]
                for (i, (m, design)) = enumerate(zip(models, designs))
                    # a = Gray.(m() .< 0.5)

                    # Images.save(joinpath(ckptpath, "optimized_design_region_$i.png"), a)
                    # Images.save(joinpath(path, "optimized_design_region_$i.png"), a)
                end
                plotsols(sols, run_probs, (path, ckptpath))

                sol = (;
                    sparam_family(S)...,
                    optimized_designs=[m() .> 0.5 for m in models], dl,
                    params=getfield.(models, :p),
                    designs,
                    design_config, path,
                    dx,
                    study,
                )

                open(joinpath(ckptpath, "solution.json"), "w") do f
                    write(f, json(cpu(sol)))
                end
                open(joinpath(path, "solution.json"), "w") do f
                    write(f, json(cpu(sol)))
                end
                println(json(sol.T, 4))
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
            println("====\n")
        end
        println("Done in $(time() - t0) .")

    end
    if length(sol.T) > 1
        println("wavelengths may have been adjusted to facilitate simulation.")
    end
    println("T-params: ")
    println(JSON.json(sol.T, 4))
    # println(sol)
    sol
end