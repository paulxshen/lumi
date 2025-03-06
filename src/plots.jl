using CairoMakie
function plotslices(a; saturation=1, path=nothing, n=1)
    fig = Figure()
    v = maximum(abs, a)
    colormap = :seismic
    colorrange = (-v, v) ./ saturation


    for i in 1:3
        for j = 1:n
            u = selectdim(a, i, round(Int, size(a, i) * j / (n + 1)))
            axis = (; aspect=size(u, 1) / size(u, 2))
            heatmap(fig[i, j], u; colormap, colorrange, axis)
        end
    end

    display(fig)
    !isnothing(path) && CairoMakie.save(path, fig)
end
function plotfield(args...; axis=(;), rulers, labels=[], kwargs...)
    fig = args[1]
    u = args[end]
    N = ndims(u)
    if N == 2
        a = u
        axis = (; axis..., aspect=(rulers[1][end] - rulers[1][1]) / (rulers[2][end] - rulers[2][1]))
        heatmap(fig, rulers..., a; axis, kwargs...)
        for (pos, text) in labels
            text!(fig, pos...; text, align=(:center, :center))
        end
    elseif N == 3
        for i in 3:-1:1
            a = selectdim(u, i, round(Int, size(u, i) / 2))
            rs = rulers[deleteat!(collect(1:3), i)]
            axis = (; axis..., aspect=(rs[1][end] - rs[1][1]) / (rs[2][end] - rs[2][1]))
            ax = fig[1, i]
            heatmap(ax, rs..., a; axis, kwargs...)
            if i == 3
                for (pos, text) in labels
                    text!(ax, pos...; text, align=(:center, :center))
                end
            end
        end
    end
end
function plotsim(prob, sol, field_names=[:Hz], geometry_names=[:ϵ]; path=nothing)
    fig = Figure()
    rulers = prob.grid.rulers
    labels = map(vcat(prob.source_instances, prob.monitor_instances)) do x
        @unpack center, label = x.tags
        (center, label)
    end
    for (i, k) = enumerate(field_names)
        a = sol.u(k)
        a /= maximum(abs, a)
        plotfield(fig[i, 1], a; rulers=rulers[k], colormap=:seismic, colorrange=(-1, 1), axis=(; title="$k"), labels)
    end
    for (i, k) = enumerate(geometry_names)
        i += length(field_names)
        a = prob.geometry(k)[1]
        a /= maximum(abs, a)
        plotfield(fig[i, 1], a; rulers=rulers.default, colormap=:grays, colorrange=(0, 1), axis=(; title="$k"), labels)
    end

    display(fig)
    if !isnothing(path)
        CairoMakie.save(joinpath(path), fig)
        println("plot saved to $path")
    end
end