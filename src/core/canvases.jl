struct Canvas
    swaps
    bbox
    lvoid
    lsolid
    symmetries
    ratio
    params
end
Base.ndims(c::Canvas) = size(c.bbox, 1)
function Canvas(swaps, bbox, lvoid, lsolid, symmetries, λ; ratio=4, params=nothing)
    symmetries = map(symmetries) do s
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
    Canvas(swaps, bbox / λ, lvoid / λ, lsolid / λ, symmetries, ratio, params)
end
struct CanvasInstance
    swaps
    start
    ratio
    _frame
    model
end
function CanvasInstance(canvas, grid, geometry; z=nothing)
    @unpack swaps, bbox, lvoid, lsolid, symmetries, ratio, params = canvas
    @unpack rulers, deltas, F = grid
    N = ndims(canvas)
    rulers = rulers.default
    start = int.(indexof.(rulers[1:N], bbox[:, 1]))
    stop = int.(indexof.(rulers[1:N], bbox[:, 2]))
    sz = stop - start
    if !isnothing(z)
        iz = round(Int, indexof.(rulers[3], z))
        push!(start, iz)
        push!(stop, iz)
    end
    dx = deltas[1][start[1]]
    ks = collect(keys(geometry))

    _frame = map(ks) do k
        k => map(geometry[k]) do a
            upsample(a[range.(start - 1, stop)...], ratio)
        end
    end |> OrderedDict
    dl = dx / ratio
    model = Blob(sz; solid_frac=1, lsolid=lsolid / dl, lvoid=lvoid / dl, symmetries, F,)
    if !isnothing(params)
        model.p .= params
    end
    start -= 1
    CanvasInstance(swaps, start, ratio, _frame, model)
end
function apply_canvas!(c, geometry, models)
    @unpack swaps, start, ratio, _frame, model = c
    ks = collect(keys(swaps))
    if :ϵ ∈ ks
        push!(ks, :invϵ)
    end
    # m = model()
    m = models[1]()
    for k in ks
        v0, v = swaps[k]
        b = Buffer(_frame[k])
        b[:] = _frame[k]
        b = place!(b, ratio + 1, m * v + (1 - m) * v0)
        a = copy(b)
        if k == :invϵ
            # a = tensorinv(a, ratio; inv=true, tensor=true)
            a = tensorinv(a, ratio)
        else
            a = [downsample(a, ratio)]
        end
        # a=tensorinv(a,ratio)

        geometry[k] = map(geometry[k], a) do g, a
            b = Buffer(g)
            b[:] = g
            b = place!(b, start, a)
            copy(b)
        end
    end
    geometry
end