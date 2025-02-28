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
function plotfield(args...; axis=(;), kwargs...)
    u = args[end]
    N = ndims(u)
    if N == 2
        a = u
        axis = (; axis..., aspect=size(a, 1) / size(a, 2))

        heatmap(args...; axis, kwargs...)
    elseif N == 3
        for i in 1:3
            a = selectdim(a, i, round(Int, round(Int, size(a, i) / 2)))
            axis = (; axis..., aspect=size(a, 1) / size(a, 2))
            fig = args[1]
            heatmap(fig[i], args[2:end-1]..., a; axis, kwargs...)
        end
    end
    fig
end
function plotsim(prob, sol, field_names=[:Hz], geometry_names=[:ϵ]; path=nothing)
    fig = Figure()
    @unpack rulers = prob.grid
    for (i, k) = enumerate(field_names)
        a = sol.u(k)
        a /= maximum(abs, a)
        plotfield(fig[i, 1], rulers[k]..., a; colormap=:seismic, colorrange=(-1, 1), axis=(; title=k))
    end
    for (i, k) = enumerate(geometry_names)
        i += length(field_names)
        a = prob.geometry(k)
        a /= maximum(abs, a)
        plotfield(fig[i, 1], rulers.default..., a; colormap=:gray, colorrange=(0, 1), axis=(; title=k))
    end

    display(fig)
    !isnothing(path) && CairoMakie.save(path, fig)
end