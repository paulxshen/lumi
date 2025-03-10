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

function supersamplemesh(a::AbstractArray{T}, ratio=4) where {T}
    N = ndims(a)
    In = LinearAlgebra.I(N)
    # spacings = _downvec.(spacings, size(a))
    a = downsample(a, ratio) do a
        p, q = @ignore_derivatives extrema(a)
        p == q & return In / a[1]

        n = imnormal(v)
        P = n * n'
        P * mean(1 ./ a) + (In - P) / mean(a)
    end
    v = map(1:N) do i
        map(1:i) do j
            map(a) do v
                v(i, j)
            end
        end
    end
    [j <= i ? v[i][j] : v[j][i] for i = 1:N, j = 1:N]
end

function supersamplemesh(meshvals::AbstractVector{<:Tuple}, rulers; tensor=false, inv=false, big=true, z=nothing)
    N = length(rulers)
    In = LinearAlgebra.I(N)
    sz = Tuple(length.(rulers) - 1)
    rulers1 = isnothing(z) ? rulers : vcat(rulers, [z])

    δ = minimum.(diff.(rulers1)) / 10
    c = map(Base.product(rulers1...)) do p
        ps = map((-δ, +δ)) do d
            Point((p + d)...)
        end
        for (i, (m, v)) = enumerate(meshvals)
            isnothing(m) && return i
            v = sideof(ps, m) .== IN
            all(v) && return i
            v[1] != v[2] && return 0
        end
    end

    a = map(CartesianIndices(sz)) do I
        I = Tuple(I)
        ms = map(Base.product(zip(I, I + 1)...)) do I
            c[I...]
        end

        i = minimum(ms)
        if i > 0 && i == maximum(ms)
            v = meshvals[i][2]
            tensor && inv && return (In) / v
            return v
        end

        start = getindex.(rulers, I)
        stop = getindex.(rulers, I + 1)
        Δ = stop - start

        if inv && tensor && big
            start = start - Δ / 2
            stop = stop + Δ / 2
            start .= max.(start, first.(rulers))
            stop .= min.(stop, last.(rulers))
            Δ = stop - start
        end

        n = 8
        δ = Δ / n
        # hits = sort(unique(ms))
        a = map(Base.product(range.(start + δ / 2, stop - δ / 2, n)...)) do p
            p = isnothing(z) ? Point(p...) : Point(p..., z)
            for (m, v) = meshvals
                if isnothing(m) || sideof(p, m) == IN
                    return v
                end
            end
        end

        if tensor && inv
            n = imnormal(a)
            P = n * n'
            P * mean(1 ./ a) + (In - P) / mean(a)
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