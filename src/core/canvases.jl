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
Base.ndims(c::CanvasInstance) = length(c.start)
function CanvasInstance(canvas, grid, geometry; z=nothing)
    @unpack swaps, bbox, lvoid, lsolid, symmetries, ratio, params = canvas
    @unpack rulers, deltas, F = grid
    N = ndims(canvas)
    rulers = rulers.default
    deltas = deltas.default
    global _d = rulers, canvas
    start = int.(indexof.(rulers[1:N], bbox[:, 1]))
    stop = int.(indexof.(rulers[1:N], bbox[:, 2]))
    sz = stop - start
    dx = deltas[1][start[1]]

    _frame = map(keys(swaps)) do k
        k => map(geometry[k]) do a
            upsample(a[range.(start - 1, stop)...], ratio)
        end[1]
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
    N = ndims(c)
    ks = collect(keys(swaps))

    # m = model()
    m = models[1]()
    for k in ks
        v0, v = swaps[k]
        b = Buffer(_frame[k])
        b[:] = _frame[k]
        b = place!(b, m * v + (1 - m) * v0, fill(ratio + 1, N))
        _a = copy(b)

        if :ϵ == k
            _ks = (k, :invϵ)
        else
            _ks = (k,)
        end

        for k = _ks
            if k == :invϵ
                # a = tensorinv(a, ratio; inv=true, tensor=true)
                a = tensorinv(_a, ratio)
            else
                a = [downsample(_a, ratio)]
            end
            # a=tensorinv(a,ratio)

            geometry[k] = map(geometry[k], a) do g, a
                b = Buffer(g)
                b[:] = g
                b = place!(b, a, start)
                copy(b)
            end
        end
    end
    geometry
end