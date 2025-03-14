abstract type AbstractMonitor end

mutable struct Monitor <: AbstractMonitor
    λmodenums
    λsmode
    λmodes

    center
    dimensions
    frame

    tags
end

function Monitor(center, dimensions, frame; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...)
    tags = OrderedDict(tags)
    tags[:center] = center
    Monitor(λmodenums, λsmode, λmodes, center, dimensions, frame, tags)
end

Base.ndims(m::Monitor) = length(m.center)

struct PlaneMonitor <: AbstractMonitor
    dims
    q

    λmodes
    # meta
    tags
    function PlaneMonitor(q, λmodes; dims, tags...)
        new(dims, q, λmodes, tags)
    end
end



wavelengths(m::Monitor) = keys(m.λmodes)

abstract type AbstractMonitorInstance end
mutable struct MonitorInstance <: AbstractMonitorInstance
    frame
    I
    box_deltas
    plane_Is
    plane_deltas
    λmodes
    _λmodes
    tags
end
@functor MonitorInstance (λmodes, _λmodes, plane_deltas)
Base.ndims(m::MonitorInstance) = length(m.center)
area(m::MonitorInstance) = m.v
wavelengths(m::MonitorInstance) = keys(m.λmodes)
Base.length(m::MonitorInstance) = 1

function MonitorInstance(m::Monitor, λ, g, ϵ, TEMP; z=nothing, mode_solutions=nothing)
    @unpack λmodes, _λmodes, box_rulers, bbox, box_deltas, I, plane_deltas, plane_Is, = _get_λmodes(m, λ, ϵ, TEMP, mode_solutions, g; z)
    @unpack F = g
    @unpack tags, frame, = m
    frame, = F.((frame,))
    MonitorInstance(frame, I, box_deltas, plane_Is, plane_deltas, λmodes, _λmodes, tags)
end

function MonitorInstance(m::PlaneMonitor, g)
    @unpack dimensions = g
    @unpack dims, q, λmodes, tags = m
    MonitorInstance(Monitor(λmodes, dimensions / 2, -dimensions / 2, dimensions / 2, getdimsperm(dims), tags), g)
end

function field(u::Map, k, m)
    @unpack I = m
    @nograd k, m, I
    a = u(k)
    getindexf(a, I[k]...)
end
Base.getindex(a::AbstractDict, k, m::MonitorInstance) = field(a, k, m)
Base.getindex(a::NamedTuple, k, m::MonitorInstance) = field(a, k, m)

