nc = 4
function transit(a::T, b) where {T}
    n = ceil(Int, log(1.2, b / a))
    r = (b / a)^(1 / n) |> T
    # a * r .^ (0:n)
    a * r .^ [zeros(T, nc)..., (1:n)...]
end
function share(v, x)
    if v[2] > v[1]
        a = v[nc+1:end]
        vcat(v[1:nc], a + x * a / sum(a))
    else
        a = v[1:end-nc]
        vcat(a + x * a / sum(a), v[end-nc+1:end])
    end
end

function partition(d, n, nres, nmax)
    nres *= (n / nmax)
    nres = max(2.1, nres)
    c = d * nres * n
        c = ceil(Int, c)
        d / c, c
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
                vr = ruler[j-1][2]
                deleteat!(ruler, i:j-1)
                insert!(ruler, i, (a, v))
                b < _b && insert!(ruler, i + 1, (b, vr))
            end
        end
    end
    deltas = map(rulers) do ruler
        map(enumerate(ruler[1:end-1])) do (i, (a, v))
            b = ruler[i+1][1]
            d = b - a

            vl = vr = dl = dr = nothing
            if i > 1
                dl = a - ruler[i-1][1]
                vl = ruler[i-1][2]
            end

            if i < length(ruler) - 1
                dr = ruler[i+2][1] - b
                vr = ruler[i+1][2]
            end

            lgood = i == 1 || v >= vl
            rgood = i == length(ruler) - 1 || v >= vr

            Δ, n = partition(d, v, nres, nmax,)
            if lgood && rgood
                fill(Δ, n)
            elseif !lgood && !rgood
                # Δ0, c0 = partition(dl, vl, nres, nmax)
                # Δ1, c1 = partition(dr, vr, nres, nmax)
                # L = transit(Δ0, Δ)
                # R = transit(Δ1, Δ)
                # s = d - (sum(L) + sum(R))
                # if s > 0
                #     Δs, ns, rm = partition(s, v, nres, nmax, true)
                #     vcat(share(L, rm / 2), fill(Δs, ns), share(reverse(R), rm / 2))
                # else
                #     Δ, n = partition(d, max(vl, vr), nres, nmax)
                #     fill(Δ, n)
                # end
            else
                if lgood
                    dlow, vlow = dl, vl
                    dhi, vhi = dr, vr
                elseif rgood
                    dlow, vlow = dr, vr 
                    dhi, vhi = dl, vl
                end
                                Δhi, = partition(dhi, vhi, nres, nmax)
                                Δs = transit(Δhi, Δ)
                                i = searchsortedfirst(cumsum(Δs ), d) - 1
                                if i >= 2nc
                                    Δs =Δs[1:i]
                                    r = d - sum(Δs )
                                    Δs =  if r >= Δs[end]
                                        m = floor(Int, r / Δs[end])
                                        vcat(Δs , fill(r / m, m))
                                    else
                                        share(Δs , r)
                                    end
                if lgood 
                    reverse(Δs)
                end
                                else
                                    Δ, n = partition(d, vhi, nres, nmax)
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