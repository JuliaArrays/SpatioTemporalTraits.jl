@testset "spatial" begin
    @testset "no units, no time" begin
        A = NamedAxisArray(reshape(1:12, 3, 4), x = 1:3, y = 1:4);
        @test_throws ArgumentError time_axis(A)
        @test !has_timedim(A)
        @test_throws ArgumentError timedim(A)
        @test ntime(A) == 1
        @test @inferred(pixel_spacing(A)) === (1,1)
        @test @inferred(spatial_directions(A)) === ((1,0),(0,1))
        @test @inferred(spatialdims(A)) === (1, 2)
        @test @inferred(spatial_order(A)) === (:x, :y)
        @test @inferred(spatial_size(A)) === (3,4)
        @test @inferred(spatial_indices(A)) === (Base.OneTo(3), Base.OneTo(4))
        assert_timedim_last(A)
    end

    @testset "units, no time" begin
        mm = u"mm"     # in real use these should be global consts
        m = u"m"
        A = NamedAxisArray(reshape(1:12, 3, 4), x = 1mm:1mm:3mm, y = 1m:2m:7m)
        @test_throws ArgumentError time_axis(A)
        @test !has_timedim(A)
        @test_throws ArgumentError timedim(A)
        @test ntime(A) == 1
        @test @inferred(pixel_spacing(A)) === (1mm,2m)
        @test @inferred(spatial_directions(A)) === ((1mm,0m),(0mm,2m))
        @test @inferred(spatialdims(A)) === (1,2)
        @test @inferred(spatial_order(A)) === (:x,:y)
        @test @inferred(spatial_size(A)) === (3,4)
        @test @inferred(spatial_indices(A)) === (Base.OneTo(3),Base.OneTo(4))
        assert_timedim_last(A)
    end

    @testset "units, time" begin
        s = u"s" # again, global const
        A = NamedAxisArray(reshape(1:12, 3, 4), x = 1:3, time = 1s:1s:4s)
        @test @inferred(times(A)) == 1s:1s:4s
        @test has_timedim(A)
        @test timedim(A) == 2
        @test ntime(A) == 4
        @test @inferred(pixel_spacing(A)) === (1,)
        @test @inferred(spatial_directions(A)) === ((1,),)
        @test @inferred(coords_spatial(A)) === (1,)
        @test @inferred(spatialorder(A)) === (:x,)
        @test @inferred(size_spatial(A)) === (3,)
        @test @inferred(spatial_indices(A)) === (Base.OneTo(3),)
        assert_timedim_last(A)
    end

    @testset "units, time first" begin
        s = u"s" # global const
        A = NamedAxisArray(reshape(1:12, 4, 3), time = 1s:1s:4s, x = 1:3)
        @test @inferred(times(A)) === 1s:1s:4s
        @test has_timedim(A)
        @test timedim(A) == 1
        @test ntimes(A) == 4
        @test @inferred(pixel_spacing(A)) === (1,)
        @test @inferred(spatial_directions(A)) === ((1,),)
        @test @inferred(coords_spatial(A)) === (2,)
        @test @inferred(spatial_order(A)) === (:x,)
        @test @inferred(spatial_size(A)) === (3,)
        @test @inferred(spatial_indices(A)) === (Base.OneTo(3),)
        @test_throws ErrorException assert_timedim_last(A)
    end
end

using ImageCore, Colors, FixedPointNumbers, ColorVectorSpace, MappedArrays
using Test
using ImageCore.ImageColors: Pixel, NumberLike, GenericImage, GenericGrayImage

@testset "Image traits" begin
    for (B, swap) in ((rand(UInt16(1):UInt16(20), 3, 5), false),
                      (rand(Gray{Float32}, 3, 5), false),
                      (rand(RGB{Float16}, 3, 5), false),
                      (bitrand(3, 5), false),
                      (rand(UInt32, 3, 5), false),
                      (view(rand(3, 2, 5), :, 1, :), false),
                      (OffsetArray(rand(3, 5), -1:1, -2:2), false),
                      (permuteddimsview(rand(5, 3), (2, 1)), true),
                      (mappedarray(identity, permuteddimsview(rand(5, 3), (2, 1))), true),
                      (colorview(RGB, zeros(3, 5), zeroarray, zeros(3, 5)), false))
        @test pixel_spacing(B) == (1,1)
        # Deprecated
        @test pixelspacing(B) == (1,1)
        if !isa(B, SubArray)
            @test spatial_directions(B) == (swap ? ((0,1),(1,0)) : ((1,0),(0,1)))
        else
            @test spatial_directions(B) == ((1,0,0), (0,0,1))
        end
        # Deprecated
        if !isa(B, SubArray)
            @test spacedirections(B) == (swap ? ((0,1),(1,0)) : ((1,0),(0,1)))
        else
            @test spacedirections(B) == ((1,0,0), (0,0,1))
        end

        @test sdims(B) == 2

        @test spatialdims(B) == (swap ? (2,1) : (1,2))
        # Deprecated
        @test coords_spatial(B) == (swap ? (2,1) : (1,2))

        @test nimages(B) == 1

        @test spatial_size(B) == (3,5)
        # Deprecated
        @test size_spatial(B) == (3,5)

        if isa(B, OffsetArray)
            @test spatial_axes(B) == (-1:1, -2:2)
        else
            @test spatial_axes(B) == (Base.OneTo(3), Base.OneTo(5))
        end
        # Deprecated
        if isa(B, OffsetArray)
            @test indices_spatial(B) == (-1:1, -2:2)
        else
            @test indices_spatial(B) == (Base.OneTo(3), Base.OneTo(5))
        end

        assert_timedim_last(B)
        @test width(B) == 5
        @test height(B) == 3
    end
end

# delibrately written in a redundant way
@testset "*Like traits" begin
    @testset "Pixel" begin
        @test NumberLike <: Pixel
        @test Number <: Pixel
        @test Gray <: Pixel
        @test RGB <: Pixel

        @test isa(oneunit(Gray), Pixel)
        @test isa(RGB(1.0, 0.0, 0.0), Pixel)
    end

    @testset "NumberLike" begin
        @test Number <: NumberLike
        @test Real <: NumberLike
        @test AbstractFloat <: NumberLike
        @test FixedPoint <: NumberLike
        @test Integer <: NumberLike
        @test Bool <: NumberLike

        @test Gray <: NumberLike
        @test Gray{<:AbstractFloat} <: NumberLike
        @test Gray{<:Bool} <: NumberLike

        @test isa(oneunit(Gray), NumberLike)
    end

    @testset "GenericImage" begin
        @test GenericGrayImage <: GenericImage
        for sz in [(3, 3), (3, 3, 3)]
            @test isa(rand(Bool, sz), GenericImage)
            @test isa(rand(N0f8, sz), GenericImage)
            @test isa(rand(Float32, sz), GenericImage)

            @test isa(rand(Gray, sz), GenericImage)
            @test isa(rand(Gray{Bool}, sz), GenericImage)
            @test isa(rand(Gray{N0f8}, sz), GenericImage)
            @test isa(rand(Gray{Float32}, sz), GenericImage)

            @test isa(rand(GrayA, sz), GenericImage)
            @test isa(rand(GrayA{N0f8}, sz), GenericImage)
            @test isa(rand(GrayA{Float32}, sz), GenericImage)

            @test isa(rand(RGB, sz), GenericImage)
            @test isa(rand(RGB{N0f8}, sz), GenericImage)
            @test isa(rand(RGB{Float32}, sz), GenericImage)

            @test isa(rand(RGBA, sz), GenericImage)
            @test isa(rand(RGBA{N0f8}, sz), GenericImage)
            @test isa(rand(RGBA{Float32}, sz), GenericImage)

            @test isa(rand(Lab, sz), GenericImage)
            @test isa(rand(Lab{Float32}, sz), GenericImage)
        end
    end

    @testset "GrayImage" begin
        for sz in [(3, 3), (3, 3, 3)]
            @test isa(rand(Bool, sz), GenericGrayImage)
            @test isa(rand(N0f8, sz), GenericGrayImage)
            @test isa(rand(Float32, sz), GenericGrayImage)

            @test isa(rand(Gray, sz), GenericGrayImage)
            @test isa(rand(Gray{Bool}, sz), GenericGrayImage)
            @test isa(rand(Gray{N0f8}, sz), GenericGrayImage)
            @test isa(rand(Gray{Float32}, sz), GenericGrayImage)
        end
    end

    @testset "dispatch" begin
        begin
            whatis(::GenericImage) = "GenericImage"
            whatis(::GenericGrayImage) = "GenericGrayImage"
            whatis(::GenericImage{<:AbstractRGB}) = "GenericRGBImage"

            whatis(::GenericImage{<:Pixel, 2}) = "Generic2dImage"
            whatis(::GenericGrayImage{<:NumberLike, 2}) = "Gray2dImage"
            whatis(::GenericImage{<:AbstractRGB, 2}) = "RGB2dImage"

            @test whatis(rand(Lab, 2, 2, 2)) == "GenericImage"

            @test whatis(rand(Lab, 2, 2)) == "Generic2dImage"

            @test whatis(rand(2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(N0f8, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Bool, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Float32, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Int64, 2, 2, 2)) == "GenericGrayImage"

            @test whatis(rand(Gray, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Gray{N0f8}, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Gray{Bool}, 2, 2, 2)) == "GenericGrayImage"
            @test whatis(rand(Gray{Float32}, 2, 2, 2)) == "GenericGrayImage"

            @test whatis(rand(2, 2)) == "Gray2dImage"
            @test whatis(rand(N0f8, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Bool, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Float32, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Int64, 2, 2)) == "Gray2dImage"

            @test whatis(rand(Gray, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Gray{N0f8}, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Gray{Bool}, 2, 2)) == "Gray2dImage"
            @test whatis(rand(Gray{Float32}, 2, 2)) == "Gray2dImage"

            @test whatis(rand(RGB, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(RGB{N0f8}, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(RGB{Float32}, 2, 2, 2)) == "GenericRGBImage"

            @test whatis(rand(BGR, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(BGR{N0f8}, 2, 2, 2)) == "GenericRGBImage"
            @test whatis(rand(BGR{Float32}, 2, 2, 2)) == "GenericRGBImage"

            @test whatis(rand(RGB, 2, 2)) == "RGB2dImage"
            @test whatis(rand(RGB{N0f8}, 2, 2)) == "RGB2dImage"
            @test whatis(rand(RGB{Float32}, 2, 2)) == "RGB2dImage"

            @test whatis(rand(BGR, 2, 2)) == "RGB2dImage"
            @test whatis(rand(BGR{N0f8}, 2, 2)) == "RGB2dImage"
            @test whatis(rand(BGR{Float32}, 2, 2)) == "RGB2dImage"
        end
    end
end

struct RowVector{T,P} <: AbstractVector{T}
    v::Vector{T}
    p::P
end

AxisIndices.has_dimnames(::Type{<:RowVector}) = true
AxisIndices.dimnames(::Type{<:RowVector}) = (:row,)
Base.axes(rv::RowVector) = axes(rv.v)


@testset "Trait Interface" begin
    img = reshape(1:24, 2,3,4)

    @test @inferred(named_axes(img)) == NamedTuple{(:dim_1, :dim_2, :dim_3)}(axes(img))
    # Deprecated
    @test @inferred(namedaxes(img)) == NamedTuple{(:dim_1, :dim_2, :dim_3)}(axes(img))

    rv = RowVector([1:10...], Dict{String,Any}())

    @test @inferred(named_axes(rv)) == NamedTuple{(:row,)}((Base.OneTo(10),))
    # Deprecated
    @test @inferred(namedaxes(rv)) == NamedTuple{(:row,)}((Base.OneTo(10),))
end
