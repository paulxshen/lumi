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

function Monitor(center, dimensions, frame, λ=1; λmodenums=nothing, λsmode=nothing, λmodes=nothing, tags...)
    center /= λ
    dimensions /= λ
    λmodenums = kmap(x -> x / λ, identity, λmodenums)
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

    center
    frame
    I
    plane_Is
    plane_deltas
    λmodes
    _λmodes
    tags
end
@functor MonitorInstance (λmodes, _λmodes, deltas)
Base.ndims(m::MonitorInstance) = length(m.center)
area(m::MonitorInstance) = m.v
wavelengths(m::MonitorInstance) = keys(m.λmodes)
Base.length(m::MonitorInstance) = 1

function MonitorInstance(m::Monitor, g, ϵ, TEMP, mode_solutions=nothing)
    @unpack λmodes, _λmodes, plane_rulers, bbox, plane_deltas, I, plane_points, plane_Is, labelpos = _get_λmodes(m, ϵ, TEMP, mode_solutions, g)
    # println("")
    MonitorInstance(m.center, m.frame, I, plane_Is, plane_deltas, λmodes, _λmodes, m.tags)
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

