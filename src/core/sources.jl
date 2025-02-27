
function completexyz(mode, N)
    # length(mode) == N && return mode
    OrderedDict([
        begin
            i = findfirst(keys(mode)) do k
                endswith(string(k), s)
            end
            F = first(string(first(keys(mode))))
            # @show keys(mode), s, i
            if isnothing(i)
                (Symbol("$F$s") => 0)
            else
                (keys(mode)[i] => mode[keys(mode)[i]])
            end
        end for s = "xyz"[1:N]
    ])
end

"""
"""
mutable struct Source
    λmodenums
    λsmode
    λmodes

    center
    dimensions
    frame
    # dimsperm

    tags
end

function Source(center, dimensions, frame, λ=1; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...)
    center /= λ
    dimensions /= λ
    λmodenums = kvmap((k, v) -> (k / λ, v), λmodenums)
    tags = OrderedDict(tags)

    Source(λmodenums, λsmode, λmodes, center, dimensions, frame, tags)
end

Base.ndims(m::Source) = length(m.center)
# isortho(m::Source) = !isnothing(m.dimsperm)

"""
    function PlaneWave(f, dims; mode...)

Constructs plane wave source

Args
- f: time function
- dims: eg -1 for wave coming from -x face
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct PlaneWave
    λmodes
    dims
    tags
    function PlaneWave(sigmodes, dims; tags...)
        new(sigmodes, dims, tags)
    end
end
@functor PlaneWave
"""
    function GaussianBeam(f, σ, c, dims; mode...)

Constructs gaussian beam source

Args
- f: time function
- σ: std dev length
- dims: eg 1 for x direction
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""
struct GaussianBeam
    f
    σ
    mode
    c
    dims
    function GaussianBeam(f, σ, c, dims; mode...)
        new(f, σ, mode, c, dims)
    end
end
@functor GaussianBeam

"""
    function Source(f, c, lb, ub, tags=""; mode...)
    function Source(f, c, dimensions, tags=""; mode...)
        
Constructs custom  source. Can be used to specify uniform or modal sources

Args
- f: time function
- c: g.lb or center of source
- lb: lower bounds wrt to c
- ub: upper bounds wrt to c
- dimensions: source dimensions in [wavelengths]
- mode: which mode to excite & their scaling constants (typically a current source, eg Jz=1)
"""

# mode(m::Source) = m.mode
# Base.string(m::Union{Source,Source}) =
#     """
#     $(m.tags): $(count((m.ub.-m.lb).!=0))-dimensional source, centered at $(m.center|>d2), spanning from $(m.lb|>d2) to $(m.ub|>d2) relative to center,"""
#  exciting $(join(keys(m.mode),", "))"""


struct SourceInstance
    sigmodes
    # o
    center
    # dimsperm
    tags
end
@functor SourceInstance (sigmodes,)
Base.ndims(m::SourceInstance) = length(m.center)

function SourceInstance(s::PlaneWave, g)
    @unpack dimensions = g
    @unpack dims, sigmodes, tags = s
    SourceInstance(Source(sigmodes, dimensions / 2, -dimensions / 2, dimensions / 2, getdimsperm(dims), tags), g)
end

function SourceInstance(s::Source, g, ϵ, TEMP, mode_solutions=nothing)
    @unpack tags, frame = s
    @unpack F, sizes = g
    C = complex(F)
    N = ndims(s)

    @unpack λmodes, _λmodes, box_size, bbox, box_deltas, I, plane_points, plane_Is, labelpos =
        _get_λmodes(s, ϵ, TEMP, mode_solutions, g)

    λs = @ignore_derivatives Array(keys(λmodes))
    modess = values(λmodes)
    if all(x -> x === (modess[1]), modess)
        iss = [eachindex(λs)]
        println("all modes are the same")
    else
        iss = cluster(λs)
    end
    sigmodes = map(iss) do is
        f = t -> sum(getindex.((Array(λs),), Array(is))) do λ
            cispi(2t / λ) |> C
        end
        _modess = getindex.((modess,), is)
        modes = _modess[round(length(_modess) / 2 + 0.1)]
        f, modes
    end

    sigmodes = reduce(vcat, [collect(zip(fill(f, length(modes)), modes)) for (f, modes) = sigmodes])

    sigmodes = [
        begin
            _f = if isa(sig, Number)
                t -> cispi(2t / sig)
            else
                sig
            end
            f = x -> convert(C, _f(x))

            mode = OrderedDict([
                begin
                    a = zeros(C, box_size)
                    for (start, v) = zip(plane_Is, v)
                        place!(a, v, start)
                    end
                    k => a
                end for (k, v) = pairs(mode)])
            mode = completexyz(mode, N)

            mode = packxyz(mode)
            mode = vmap(mode) do v
                frame * v
            end
            mode = unpackxyz(mode)
            global _gf2 = mode, I
            mode = namedtuple([
                k => begin
                    v = mode(k)
                    @assert !any(isnan, v)

                    if all(iszero, v)
                        a = 0
                    else
                        a = zeros(C, sizes[k])
                        place!(a, v, first.(I[k]))
                    end

                end for k = sort(keys(mode))])
            global _gf3 = mode
            (f, mode)
        end for (sig, mode) = sigmodes]
    # @show center, g.lb, labelpos

    SourceInstance(sigmodes, labelpos, tags)
end

# Complex
function (s::SourceInstance)(t::Real)
    namedtuple([
        k => sum(s.sigmodes) do (f, _mode)
            real(f(t) .* _mode[k])
        end for k = keys(s.sigmodes[1][2])])
end

function E2J(d)
    k0 = filter(k -> string(k)[1] == 'E', keys(d))
    k = [Symbol("J" * string(k)[2:end]) for k in k0]
    v = [d[k] for k in k0]
    dict(Pair.(k, v))
end
function EH2JM(d::T) where {T}
    dict(Pair.(replace(keys(d), :Ex => :Jx, :Ey => :Jy, :Ez => :Jz, :Hx => :Mx, :Hy => :My, :Hz => :Mz), values(d)))
end

function _get_λmodes(sm, ϵ, TEMP, mode_solutions, g)
    @unpack center, dimensions, tags, frame, λmodenums, λmodes, λsmode = sm
    @unpack F, sizes, rulers, all_field_names, offsets, dx = g
    N = ndims(sm)

    C = complex(g.F)

    if isa(sm, Source)
        ks = [k for k = all_field_names if string(k)[1] ∈ ('J')]
    elseif isa(sm, Monitor)
        ks = [k for k = all_field_names if string(k)[1] ∈ ('E', 'H')]
    end

    v = frame * [dimensions..., 0]
    rulers = rulers[:default]
    start = indexof.(rulers, center - v / 2)
    stop = indexof.(rulers, center + v / 2)

    s = Int.(sign.(v))
    s += s .== 0
    start = s .* (floor.(Int, s .* start))
    stop = s .* (ceil.(Int, s .* stop))
    box_size = Tuple(abs.(stop - start))

    s = v .>= 0
    corner = s + .!s .* box_size
    signed_bbox = [getindex.(rulers, start) getindex.(rulers, stop)]
    signed_plane_rulers = getindex.(rulers, map(start, stop) do a, b
        a == b ? (a:b) : a:sign(b - a):b
    end)
    signed_plane_rulers -= first.(signed_plane_rulers)

    start, stop = min.(start, stop), max.(start, stop)
    bbox = [getindex.(rulers, start) getindex.(rulers, stop)]
    box_rulers = getindex.(rulers, range.(start, stop))
    box_deltas = diff.(box_rulers)
    plane_rulers = box_rulers - first.(box_rulers)

    global I = dict([k => begin
        l, r = eachcol(offsets[k])
        range.(start - l, stop - r - 1, box_size)
    end for k = ks])
    labelpos = round(Int, indexof.(rulers, center))

    # 
    if !isnothing(λmodenums)
        global plane_size = floor.(Int, dimensions / dx)
        global P = dx * frame[:, 1:end-1]
        global signed_plane_start = center - P * collect((plane_size - 1) / 2)

        global plane_points = map(CartesianIndices(Tuple(plane_size))) do I
            P * collect(Tuple(I) - 1) + signed_plane_start
        end
        global plane_Is = map(plane_points) do p
            v = indexof.(rulers, p) - start + F(0.5)
            v = max.(v, 1)
            v = min.(v, box_size)
        end
        plane_deltas = dx

        global ϵmode = samplemesh(ϵ, plane_points) .|> F
        λmodes = OrderedDict([λ => begin
            modes = solvemodes(ϵmode, dx, λ, maximum(mns) + 1, TEMP; mode_solutions)[mns+1]
            map(modes) do mode
                if isa(sm, Monitor)
                    ks = filter(k -> string(k)[1] ∈ ('E', 'H'), keys(mode))
                else
                    ks = filter(k -> string(k)[1] ∈ ('J'), keys(mode))
                end
                namedtuple(ks .=> mode.(ks))
            end
        end for (λ, mns) = pairs(λmodenums)])
        # display(heatmap(ϵmode))

        if isa(sm, Monitor)
            λmodes = vmap(λmodes) do modes
                normalize_mode.(modes, (plane_deltas,))
            end
            _λmodes = vmap(λmodes) do modes
                mirror_mode.(modes)
            end
        else
            _λmodes = nothing
        end
    elseif !isnothing(λsmode)
        # λs, mode = λsmode
        # λs = F.(λs)
        # mode = vmap(Symbol, identity, mode)
        # mode = OrderedDict([
        #     begin
        #         v = mode[k]
        #         if isa(v, Number)
        #             v = v * block
        #         end
        #         k => v
        #     end for k = keys(mode) if k ∈ ks])

        # if isa(sm, Monitor)
        #     mode = normalize_mode(mode, md)
        #     _mode = mirror_mode(mode; flip=false)
        # end

        # λmodes = OrderedDict([λ => [mode] for λ = λs])
        # if isa(sm, Monitor)
        #     _λmodes = OrderedDict([λ => [_mode] for λ = λs])
        # else
        #     _λmodes = nothing
        # end
    end

    λmodes = fmap(F, λmodes)
    λmodes = sort(λmodes, by=kv -> kv[1])

    if !isnothing(_λmodes)

        _λmodes = fmap(F, _λmodes)
        _λmodes = sort(_λmodes, by=kv -> kv[1])
    end
    # m = vmap(x -> C.(x), m)

    global r = (; λmodes, _λmodes, box_size, box_rulers, bbox, box_deltas, I, plane_points, plane_Is, plane_deltas, labelpos)
end