function transit(a::T, b) where {T}
    n = ceil(Int, log(1.4, b / a))
    r = (b / a)^(1 / n) |> T
    # a * r .^ (0:n)
    a * r .^ [0, 0, (1:n)...]
end
function share(v, x)
    if v[2] > v[1]
        a = v[2:end]
        vcat([v[1]], a + x * a / sum(a))
    else
        a = v[1:end-1]
        vcat(a + x * a / sum(a), [v[end]])
    end
end

function partition(d, n, nres, nmax, withrem=false)
    nres *= (n / nmax)
    nres = max(2.1, nres)
    c = d * nres * n
    if withrem
        c0 = c
        c = floor(Int, c)
        r = c0 - c
        d / c0, c, r * d / c0
    else
        c = ceil(Int, c)
        d / c, c
    end
end
function makemesh(mvs, bbox, nres)
    nmax = maximum(last.(mvs))

    rulers = [[(a, last(mvs)[2]), (b, nothing)] for (a, b) = eachrow(bbox)]
    for (m, v) = reverse(mvs[1:end-1])
        x = boundingbox(m)
        @debug x
        a = ustrip.(getfield(coords(x.min), :coords))
        b = ustrip.(getfield(coords(x.max), :coords))
        for (ruler, a, b) = zip(rulers, a, b)
            i = searchsortedfirst(ruler, a; by=first)
            j = searchsortedfirst(ruler, b; by=first)
            if i <= length(ruler)
                j = min(j, length(ruler))
                _b = ruler[j][1]
                v1 = ruler[j-1][2]
                deleteat!(ruler, i:j-1)
                insert!(ruler, i, (a, v))
                b < _b && insert!(ruler, i + 1, (b, v1))
            end
        end
    end
    deltas = map(rulers) do ruler
        map(enumerate(ruler[1:end-1])) do (i, (a, v))
            b = ruler[i+1][1]
            d = b - a

            v0 = v1 = d0 = d1 = nothing
            if i > 1
                d0 = a - ruler[i-1][1]
                v0 = ruler[i-1][2]
            end

            if i < length(ruler) - 1
                d1 = ruler[i+2][1] - b
                v1 = ruler[i+1][2]
            end

            l = i == 1 || v >= v0
            r = i == length(ruler) - 1 || v >= v1

            Δ, n = partition(d, v, nres, nmax,)
            if l && r
                fill(Δ, n)
            elseif !l && !r
                Δ0, c0 = partition(d0, v0, nres, nmax)
                Δ1, c1 = partition(d1, v1, nres, nmax)
                L = transit(Δ0, Δ)
                R = transit(Δ1, Δ)
                s = d - (sum(L) + sum(R))
                if s > 0
                    Δs, ns, rm = partition(s, v, nres, nmax, true)
                    vcat(share(L, rm / 2), fill(Δs, ns), share(reverse(R), rm / 2))
                else
                    Δ, n = partition(d, max(v0, v1), nres, nmax)
                    fill(Δ, n)
                end
            elseif l
                Δ1, c1 = partition(d1, v1, nres, nmax)
                R = transit(Δ1, Δ)
                i = searchsortedfirst(cumsum(R), d) - 1
                if i >= 3
                    R = R[1:i]
                    s = d - sum(R)
                    Δs, ns, rm = partition(s, v, nres, nmax, true)
                    vcat(fill(Δs, ns), share(reverse(R), rm))
                else
                    Δ, n = partition(d, v1, nres, nmax)
                    fill(Δ, n)
                end
            else
                Δ0, c0 = partition(d0, v0, nres, nmax)
                L = transit(Δ0, Δ)
                i = searchsortedfirst(cumsum(L), d) - 1
                if i >= 3
                    L = L[1:i]
                    s = d - sum(L)
                    Δs, ns, rm = partition(s, v, nres, nmax, true)
                    vcat(share(L, rm), fill(Δs, ns))
                else
                    Δ, n = partition(d, v0, nres, nmax)
                    fill(Δ, n)
                end
            end
        end
    end
    r = map(zip(rulers, deltas)) do (r, v)
        start = r[1][1]
        [start, (start + cumsum(reduce(vcat, v)))...]
    end
end

function samplemesh(meshvals, points; z=nothing)
    map(points) do point
        if !isnothing(z)
            point = Point(point..., z)
        else
            point = Point(point...)
        end
        for (i, (m, v)) = enumerate(meshvals)
            if isnothing(m) || sideof(point, m) == IN
                return v
            end
        end
    end
end