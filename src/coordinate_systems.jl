
"""
    CoordinateSystem

A trait for representing coordinate systems.
"""
abstract type CoordinateSystem end

struct ImageCoordinates <: CoordinateSystem end

@inline function CoordinateSystem(::Type{T}) where {T}
    if parent_type(T) <: T
        return ImageCoordinates()
    else
        return CoordinateSystem(parent_type(T))
    end
end
