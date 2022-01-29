
"""
    CoordinateSystem(:Type{T})

* Required methods
    * `Base.axes(::CoordinateSystem)`
    * `ArrayInterface.axes_types(::Type{CoordinateSystem})`
* Optional methods
    * `Base.axes(::CoordinateSystem, dim::Int)`
    * `Base.permutedims(::CoordinateSystem{1})`
    * `Base.permutedims(::CoordinateSystem{2})`
    * `Base.permutedims(::CoordinateSystem{N}, dims::Tuple{Vararg{Any,N}})`
"""
abstract type CoordinateSystem{N} end

Base.ndims(::Type{<:CoordinateSystem{N}}) where {N} = N::Int

# TODO document that trailing axes are Cartesian unless otherwise specified
@inline Base.axes(x::CoordinateSystem, dim::StaticInt{D}) where {D} = axes(x, D)
@inline function Base.axes(x::CoordinateSystem{N}, dim::Int) where {N}
    if dim > N
        return CartesianAxis()
    else
        return getfield(axes(x), dim)
    end
end

@inline function Base.permutedims(x::CoordinateSystem{1})
    CartesianSystem(CartesianAxis(), getfield(axes(x), 1))
end
@inline Base.permutedims(x::CoordinateSystem{2}) = CartesianSystem(reverse(axes(x)))
@inline function Base.permutedims(x::CoordinateSystem{N}, dims::Tuple{Vararg{Any,N}}) where {N}
    CartesianSystem(map(d -> axes(x, d), dims))
end

"""
    EulerAxis

Represents an X, Y, or Z euler axis.
"""
struct EulerAxis{N} <: CoordinateSystem{1}
    name::N

    EulerAxis{StaticSymbol{:X}}(::StaticSymbol{:X}) = new{StaticSymbol{:X}}(StaticSymbol{:X}())
    EulerAxis{StaticSymbol{:Y}}(::StaticSymbol{:Y}) = new{StaticSymbol{:Y}}(StaticSymbol{:Y}())
    EulerAxis{StaticSymbol{:Z}}(::StaticSymbol{:Z}) = new{StaticSymbol{:Z}}(StaticSymbol{:Z}())
    EulerAxis{StaticSymbol{:X}}() = EulerAxis{StaticSymbol{:X}}(StaticSymbol{:X}())
    EulerAxis{StaticSymbol{:Y}}() = EulerAxis{StaticSymbol{:Y}}(StaticSymbol{:Y}())
    EulerAxis{StaticSymbol{:Z}}() = EulerAxis{StaticSymbol{:Z}}(StaticSymbol{:Z}())
    function EulerAxis{Symbol}(s::Symbol)
        (s === :X || s === :Y || s === :Z) && return new{Symbol}(s)
        error("EulerAxis must be :X, :Y, or :Z")
    end
    EulerAxis(x::Union{Symbol,StaticSymbol}) = EulerAxis{typeof(x)}(x)
end
const XAxis = EulerAxis{StaticSymbol{:X}}
const YAxis = EulerAxis{StaticSymbol{:Y}}
const ZAxis = EulerAxis{StaticSymbol{:Z}}

Base.axes(x::EulerAxis) = (x,)

"""
    TemporalAxis

Represents a temporal axis.
"""
struct TemporalAxis <: CoordinateSystem{1} end

Base.axes(x::TemporalAxis) = (x,)

"""
    CartesianAxis

Represents a generic axis of a cartesian coordinate system.
"""
struct CartesianAxis <: CoordinateSystem{1} end

Base.axes(x::CartesianAxis) = (x,)

"""
    CartesianSystem
    
A subtype of `CoordinateSystem` composed of 
"""
struct CartesianSystem{N,A<:Tuple{Vararg{<:CoordinateSystem{1},N}}} <: CoordinateSystem{N}
    axes::A
end

Base.axes(x::CartesianSystem) = getfield(x, :axes)

## Arrays -> CoordinateSystem
@inline function CoordinateSystem(::Type{T}) where {T}
    if parent_type(T) <: T
        return CartesianSystem(_maybe_eueler(map(CoordinateSystem, dimnames(T))))
    else
        return CoordinateSystem(parent_type(T))
    end
end
_maybe_euler(x) = x
_maybe_euler(::Tuple{}) = CartesianAxis()
_maybe_euler(::NTuple{1,StaticSymbol{:_}}) = CartesianSystem((XAxis(),))
_maybe_euler(::NTuple{2,StaticSymbol{:_}}) = CartesianSystem((XAxis(),YAxis()))
_maybe_euler(::NTuple{3,StaticSymbol{:_}}) = CartesianSystem((XAxis(),YAxis(), ZAxis()))
function _maybe_euler(::NTuple{N,StaticSymbol{:_}}) where {N}
    CartesianSystem((XAxis(),YAxis(), ZAxis(), ntuple(CartesianAxis(), Val(N - 3))))
end

CoordinateSystem(x::Union{Adjoint,Transpose}) = permutedims(CoordinateSystem(parent(x)))
@inline function CoordinateSystem(x::PermutedDimsArray{T,N,I}) where {T,N,I}
    permutdims(CoordinateSystem(parent(x)), static(I))
end
@inline function CoordinateSystem(x::SubArray)
    CartesianSystem(Static.permute(axes(CoordinateSystem(parent(x))), to_parent_dims(x)))
end
if isdefined(Base, :ReshapedReinterpretArray)
    @inline function CoordinateSystem(A::Base.ReshapedReinterpretArray{T,N,S}) where {T,N,S}
        if sizeof(S) > sizeof(T)
            return CartesianSystem(CoordinateSystem(T), axes(CoordinateSystem(parent(A)))...)
        elseif sizeof(S) < sizeof(T)
            return tail(axes(CoordinateSystem(parent(A))))
        else
            return CoordinateSystem(parent(A))
        end
    end
end

CoordinateSystem(::Type{Union{StaticSymbol{:Time},StaticSymbol{:time}}}) = TemporalAxis()
CoordinateSystem(::Type{<:Dates.AbstractTime}) = TemporalAxis()
