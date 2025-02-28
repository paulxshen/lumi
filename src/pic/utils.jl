

# function lumicu()::Cint
#     if !isempty(ARGS)
#         picrun(ARGS[1]; array=cu)
#     end
#     return 0
# end

function lastrun(name=nothing; study=nothing, wd=joinpath(pwd(), "runs"))
    if !isnothing(name)
        p = joinpath(wd, name)
        println("using simulation folder $p")
        return p
    end

    l = filter(isdir, readdir(wd, join=true))
    sort!(l, by=p -> Dates.unix2datetime(mtime(p)), rev=true)

    if !isnothing(study)
        for p = l
            try
                if open(joinpath(p, "solution.json")) do f
                    JSON.parse(f)["study"]
                end == study
                    println("using simulation folder $p")
                    return p
                end
            catch e
                println(e)
            end
        end
    end
    p = l[1]
    println("using simulation folder $p")
    return p
end

function calc_sparams(run_probs,)
    lminloss = 0
    sols = solve.(run_probs)

    ignore_derivatives() do
    end

    # return sols[1]
    coeffs = OrderedDict()
    for (sol, prob) in zip(sols, run_probs)
        @unpack source_instances, monitor_instances = prob
        source_port = source_instances(1).tags.port
        source_mn = first(source_instances).tags.wavelength_mode_numbers(1)[1]
        for (i, m) = enumerate(monitor_instances)
            for (w, λ) = enumerate(keys(m.tags.wavelength_mode_numbers))
                # _λ = @ignore_derivatives Base.round(parse(F, λ), digits=4)
                for mn = m.tags.wavelength_mode_numbers[λ]
                    monitor_port = m.tags.port
                    if !haskey(coeffs, λ)
                        coeffs[λ] = OrderedDict()
                    end
                    s = "$monitor_port@$mn," * "$source_port@$source_mn"
                    s = Symbol(s)
                    coeffs[λ][s] = (sol("a+", i, w, mn), sol("a-", i, w, mn))
                end
            end
        end
    end
    # return coeffs(1)(1)[1] |> abs2

    S = OrderedDict([λ => OrderedDict([k => begin
        s = ignore() do
            split(string(k), ",")[2]
        end
        # Symbol(
        coeffs[λ][k][1] / coeffs[λ][Symbol("$s,$s")][2]
    end for (k) = keys(coeffs[λ])]) for (λ) = keys(coeffs)])
    # if source_mn == mn == 0
    #     coeffs[λ]["$monitor_port,$source_port")] = v
    # end
    return (; S, sols, lminloss)
end
function make_geometry(masks, margins, lb, dl, geometry, designs, design_config, materials; ϵeff=nothing, F=Float32, perturb=nothing)
    isnothing(masks) && return geometry
    mat = design_config.fill.material
    namedtuple([k => begin
        a = geometry[k]
        T = typeof(a)
        if k == :ϵ
            f = @ignore_derivatives materials(mat)(:epsilon) |> F
            v = @ignore_derivatives minimum(a)
            if perturb == k
                f *= convert.(F, 1.001)
            end

            b = Zygote.Buffer(a)
            copyto!(b, a)

            for (mask, margin, design) in zip(masks, margins, designs)

                o = round.(Int, (design.bbox[1] - lb) / dl) + 1
                if ndims(b) == 3
                    o = [o..., 1 + round(Int, (zcore - zmin) / dl)]
                    mask = stack(fill(mask, round(Int, thickness / dl)))
                else
                    # f = ϵeff[2] |> F
                end
                o -= margin
                mask = mask * (f - v) + v
                # b[range.(o, o .+ size(mask) .- 1)...] = mask
                b[[i:j for (i, j) = zip(o, o .+ size(mask) .- 1)]...] = mask
            end
            copy(b) |> constructor(T)
        else
            # println("no fill for $k")
            a
        end
    end for k = keys(geometry)])
end
# using GLMakie: volume
function vis(sol, prob, field=:Hz, path=".")
    @unpack u, p, _p = sol |> cpu
    prob = prob |> cpu
    @unpack monitor_instances, source_instances, λ, = prob
    @unpack deltas, spacings, bbox, dl = prob.grid
    u = u(field)
    N = ndims(u)
    # if N == 3
    #     volume(u) |> display
    #     heatmap(u[:, :, round(Int, size(u, 3) / 2)]) |> display
    # else
    #     heatmap(u) |> display
    # end
    # return
    g = _p.ϵ
    u = upsample(u, spacings)
    # g = imresize(g, size(u))
    bbox /= dl
    ratio = int(deltas[1] / dl)

    plt = quickie(u, g; dl, λ, monitor_instances, ratio, source_instances, bbox)
    display(plt)

    try
        CairoMakie.save(joinpath(path, "run.png"), plt,)
    catch e
        println("save plot failed")
        println(e)
    end
end
function plotsols(sols, probs, path,)
    # try
    for (i, (prob, sol)) in enumerate(zip(probs, sols))
        # try
        @unpack u, p = sol |> cpu
        prob = prob |> cpu
        a = u(:Hz)
        heatmap(a[:, :, round(Int, size(a, 3) / 2)]) |> display
        #     @unpack monitor_instances, source_instances, λ, = prob
        #     @unpack deltas, spacings, bbox, dl = prob.grid
        #     u = u.H.Hz
        #     N = ndims(u)
        #     # if N == 3
        #     #     volume(u) |> display
        #     #     heatmap(u[:, :, round(Int, size(u, 3) / 2)]) |> display
        #     # else
        #     #     heatmap(u) |> display
        #     # end
        #     # return
        #     g = _p.ϵ
        #     u = upsample(u, spacings)
        #     # g = imresize(g, size(u))
        #     bbox /= dl
        #     ratio = int(deltas[1] / dl)

        #     plt = quickie(u, g; dl, λ, monitor_instances, ratio, source_instances, bbox)
        #     display(plt)

        #     if !isa(path, Base.AbstractVecOrTuple)
        #         path = (path,)
        #     end
        #     for path = path
        #         try
        #             CairoMakie.save(joinpath(path, "run_$i.png"), plt,)
        #         catch e
        #             println("save plot failed")
        #             println(e)
        #         end
        #     end
        # end
        # catch e
        #     println("plot failed")
        #     println(e)
    end
end
