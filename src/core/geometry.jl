function apply_subpixel_averaging(geometry::Number, args...)
    geometry
end
function apply_subpixel_averaging(geometry, field_lims)
    namedtuple([k => begin
        a = geometry[k]
        if isa(a, Number)
            a
        else
            E = r"E.*"
            H = r"H.*"
            f = (; ϵ=E, σ=E, μ=H, m=H, invϵ=E)[k]
            lrs = field_lims(f)
            @assert length(lrs) > 0
            map(lrs) do lr
                crop(a, lr[:, 1] - 0.5, size(a) - 0.5 - lr[:, 2])
            end
        end
        # end for (k, a) = pairs(geometry)])
    end for k = keys(geometry)])
end

pad_geometry(geometry::Number, args...) = geometry
function pad_geometry(geometry, geometry_padvals, geometry_padamts)
    namedtuple([
        k => begin
            a = geometry[k]
            if a isa AbstractArray && k in keys(geometry_padvals) && k in keys(geometry_padamts)
                map(a) do a
                    pad(a, geometry_padvals[k][1], eachcol(geometry_padamts[k])...)
                end
            else
                a
            end
        end for k = keys(geometry)])
end

_size(s::Real, n) = int(s * n)
_size(s, _) = int(sum(s))

function tensorinv(a::AbstractArray{T}, field_lims, spacings) where {T}
    N = ndims(a)
    # lims = values(field_lims(r"E.*"))
    lims = @ignore_derivatives [field_lims("E$v") for v = collect("xyz")[1:N]]
    spacings = _downvec.(spacings, size(a))

    border = 0
    margin = map(spacings) do s
        (1 + border) * max(first(s), last(s))
    end
    ap = pad(a, :replicate, margin)

    v = map(1:N) do i
        map(1:i) do j
            li, ri = eachcol(lims[i])
            lj, rj = eachcol(lims[j])
            l = min.(li, lj) - 0.5
            _l = max.(li, lj) + 0.5
            Δ = _l - l
            start = round((l + 1) * margin) + 1

            # global _d = Δ, spacings, lims

            api = ap[range.(start, start + sum.(spacings) + int((Δ - 1) * last.(spacings)) - 1)...]
            ranges = [[
                range(cum - space + 1, int(Δi .* space) + cum - space)
                for (cum, space) = zip(cumsum(spacing), spacing)
            ] for (Δi, spacing) = zip(Δ, spacings)]

            downsample_by_range(api, ranges) do a
                p, q = extrema(a)
                p == q & return (i == j) / p

                n = imnormal(v)
                Pij = n[i] * n[j]
                Pij * mean(1 ./ a) + ((i == j) - Pij) / mean(a)
            end |> T
        end
    end

    [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
end

function tensorinv(meshvals::AbstractVector{<:Tuple}, rulers; tensor=false, inv=false, z=nothing)
    N = length(rulers)
    In = LinearAlgebra.I(N)
    ns = length.(rulers)
    default = meshvals[end][2]
    F = eltype(rulers[1])
    a = map(Base.product(Base.OneTo.(ns - 1)...)) do I
        start = getindex.(rulers, I)
        stop = getindex.(rulers, I .+ 1)
        Δ = stop - start

        start -= Δ / 4
        stop += Δ / 4
        start .= max.(start, first.(rulers) + F(0.001))
        stop .= min.(stop, last.(rulers) - F(0.001))
        center = (start + stop) / 2

        if !isnothing(z)
            push!(center, z)
            push!(start, z)
            push!(stop, z)
            push!(Δ, 0)
        end

        # box = Box(Point(start...), Point(stop...))
        point = Point(center...)
        hits = fill(false, length(meshvals))
        for (i, (m, v)) = enumerate(meshvals)
            # if !isnothing(m) && intersects(box, m)
            if !isnothing(m) && xor(sideof(Point(start...), m) != OUT, sideof(Point(stop...), m) != OUT)
                hits[i] = true
            elseif (!any(hits) && isnothing(m)) || (!isnothing(m) && sideof(point, m) != OUT)
                tensor && inv && return (In) / v
                return v
            end
        end

        n = 4
        δ = Δ / n
        a = map(Base.product(range.(start + δ / 2, stop - δ / 2, n)...)) do p
            for (m, v) = meshvals[hits]
                sideof(Point(p...), m) != OUT && return v
            end
            default
        end

        if tensor && inv
            n = imnormal(a)
            P = n * n'
            P * mean(1 ./ a) + (LinearAlgebra.I - P) / mean(a)
        elseif !tensor
            v = mean(a)
            inv && return 1 / v
            v
        end
    end

    if tensor && inv
        v = map(1:N) do i
            map(1:i) do j
                map(a) do v
                    v(i, j)
                end
            end
        end
        [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
    elseif !tensor && !inv
        [a]
    end
end

function downsamplefield(a::AbstractArray{T}, field_lims, spacings) where {T}
    N = ndims(a)
    lims = values(field_lims(r"E.*"))
    spacings = _downvec.(spacings, size(a))
    rulers = cumsum.(spacings)

    map(lims) do lims
        I = map(eachrow(lims), rulers) do (l, r), ruler
            map(range(l, r, length=int(r - l + 1))) do i
                start = 1 + getindexs(ruler, i) |> int
                stop = getindexs(ruler, i + 1) |> int

                start:stop
            end
        end
        downsample_by_range(mean, a, I)
    end
end