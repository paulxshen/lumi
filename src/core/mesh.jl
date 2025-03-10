nc = 4

function transit(a::T, b) where {T}
    n = ceil(Int, log(1.2, b / a))
    r = (b / a)^(1 / n) |> T
    # a * r .^ (0:n)
    a * r .^ [zeros(T, nc)..., (1:n)...]
end

function share(v, x)
    if v[end] > v[1]
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

    rulers = [[[a, last(mvs)[2]], [b, nothing]] for (a, b) = eachrow(bbox)]
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
                insert!(ruler, i, [a, v])
                b < _b && insert!(ruler, i + 1, [b, vr])
            end
        end
    end
    for v=rulers
        i=1
        while i<length(v)
            d=v[i+1][1]-v[i][1]
            if d<1f-4
                if isnothing(v[i][2])
                    v[i][2]=v[i+1][2]
                    deleteat!(v,i+1)
                else
                    deleteat!(v,i)
                end
            else
                i+=1
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
                Δl, = partition(dl, vl, nres, nmax)
                Δr,  = partition(dr, vr, nres, nmax)
                L = transit(Δl, Δ)
                R = transit(Δr, Δ)
                r= d - (sum(L) + sum(R))
                if r > 0
                    Δ= L[end]+R[end]
                    if r >=Δ
                        m = floor(Int, r / Δ)
                        Δs =       vcat(L, fill(r / m, m),reverse(R))
                    else
                        vcat(share(L, r / 2),  share(reverse(R), r / 2))
                    end
                else
                    Δ, n = partition(d, max(vl, vr), nres, nmax)
                    fill(Δ, n)
                end
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
                if i >= nc+2
                    Δs =Δs[1:i]
                    r = d - sum(Δs )
                   if r >= Δs[end]
                        m = floor(Int, r / Δs[end])
                        Δs =       vcat(Δs , fill(r / m, m))
                    else
                        Δs =   share(Δs , r)
                    end
                    if lgood 
                        Δs =   reverse(Δs)
                        else
                        Δs
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