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
function Canvas(swaps, bbox, lvoid, lsolid, symmetries; ratio=4, params=nothing)
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
    Canvas(swaps, bbox, lvoid, lsolid, symmetries, ratio, params)
end
struct CanvasInstance
    swaps
    start
    ratio
    model
end
Base.ndims(c::CanvasInstance) = length(c.start)
function CanvasInstance(canvas, λ, grid; z=nothing)
    @unpack swaps, bbox, lvoid, lsolid, symmetries, ratio, params = canvas
    @unpack rulers, deltas, F = grid

    bbox /= λ
    lvoid /= λ
    lsolid /= λ

    N = ndims(canvas)
    rulers = rulers.default
    deltas = deltas.default
    start = int.(indexof.(rulers[1:N], bbox[:, 1]), 0.1)
    stop = int.(indexof.(rulers[1:N], bbox[:, 2]), 0.1)
    sz = stop - start
    dx = deltas[1][start[1]]

    model = Blob(sz * ratio; solid_frac=1, lsolid=ratio * lsolid / dx, lvoid=ratio * lvoid / dx, symmetries, F,)
    if !isnothing(params)
        model.p .= params
    end
    CanvasInstance(swaps, start, ratio, model)
end

function apply_canvas!(c, geometry, models)
    @unpack swaps, start, ratio, model = c
    N = ndims(c)
    ks = collect(keys(swaps))

    # m = model()
    m = models[1]()
    for k in ks
        v0, v = swaps[k]
        _a = m * v + (1 - m) * v0

        if :ϵ == k
            _ks = (k, :invϵ)
        else
            _ks = (k,)
        end

        for k = _ks
            if k == :invϵ
                # a = supersamplemesh(a, ratio; inv=true, tensor=true)
                a = supersamplemesh(_a, ratio)
            else
                a = [downsample(_a, ratio)]
            end
            # a=supersamplemesh(a,ratio)
            geometry[k] = map(geometry[k], a) do g, a
                b = Buffer(g)
                b[:] = g
                b = place!(b, a, start; additive=false)
                copy(b)
            end
        end
    end
    geometry
end