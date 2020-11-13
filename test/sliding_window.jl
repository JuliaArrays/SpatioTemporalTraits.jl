
using SpatioTemporalTraits: sliding_window

@testset "SlidingWindow" begin
    axis = SimpleAxis(20)
    @test @inferred(collect(sliding_window(axis; window_size=3))) == [1:3, 4:6, 7:9, 10:12, 13:15, 16:18]
    @test @inferred(collect(sliding_window(axis; window_size=3, first_pad=1))) == [2:4, 5:7, 8:10, 11:13, 14:16, 17:19]
    @test @inferred(collect(sliding_window(axis; window_size=3, first_pad=1, last_pad=1))) == [2:4, 5:7, 8:10, 11:13, 14:16, 17:19]
    @test @inferred(collect(sliding_window(axis; window_size=3, first_pad=1, last_pad=2))) == [2:4, 5:7, 8:10, 11:13, 14:16]
    @test @inferred(collect(sliding_window(axis; window_size=3, first_pad=1, last_pad=2, dilation=2))) == [2:2:4, 5:2:7, 8:2:10, 11:2:13, 14:2:16]
    @test @inferred(collect(sliding_window(axis; window_size=3, first_pad=1, last_pad=2, stride=2))) == [2:4, 7:9, 12:14]
    axis = Axis(range(2.0, step=3.0, length=20))
    itr = sliding_window(axis; window_size=9.0)
    @test length(itr) == 6
    @test iterate(itr) == (1:3, 0)
    @test first(itr) == 1:3
    @test iterate(itr, 0) == (4:6, 3)
    @test last(itr) == 16:18
    @test @inferred(collect(itr)) == [1:3, 4:6, 7:9, 10:12, 13:15, 16:18]
    @test @inferred(collect(sliding_window(axis; window_size=9.0, first_pad=3.0))) == [2:4, 5:7, 8:10, 11:13, 14:16, 17:19]
    @test @inferred(collect(SlidingWindow(axis; window_size=9.0, first_pad=3.0, last_pad=3.0))) == [2:4, 5:7, 8:10, 11:13, 14:16, 17:19]
    @test @inferred(collect(SlidingWindow(axis; window_size=9.0, first_pad=3.0, last_pad=6.0))) == [2:4, 5:7, 8:10, 11:13, 14:16]
    itr = sliding_window(axis; window_size=9.0, first_pad=3.0, last_pad=6.0, dilation=6.0)
    @test iterate(itr) == (2:2:4, 1)
    @test iterate(itr, 1) == (5:2:7, 4)
    @test @inferred(collect(itr)) == [2:2:4, 5:2:7, 8:2:10, 11:2:13, 14:2:16]
    itr = sliding_window(axis; window_size=9.0, first_pad=3.0, last_pad=6.0, stride=6.0)
    @test iterate(itr) == (2:4, 1)
    @test iterate(itr, 1) ==(7:9, 6)
    @test @inferred(collect(itr)) == [2:4, 7:9, 12:14]
    @test iterate(itr, 11) === nothing
end

@testset "NDSlidingWindow" begin
    axis = SimpleAxis(20)
    A = CartesianIndices((axis, axis, axis));
    axsitr = sliding_window(A; window_size=(3,3,3), first_pad=(1,1,1), last_pad=(2,2,2), stride=(2,2,2))
    @test SpatioTemporalTraits.firstinc(axsitr.sliding_windows) == ((2:4, 2:4, 2:4), (1, 1, 1))
    @test iterate(axsitr) == (Iterators.ProductIterator((2:4, 2:4, 2:4)), ((2:4, 2:4, 2:4), (1, 1, 1)))
    @test first(axsitr) == Iterators.ProductIterator((2:4, 2:4, 2:4))
    @test last(axsitr) == Iterators.ProductIterator((12:14, 12:14, 12:14))
    @test iterate(axsitr, ((12:14, 12:14, 12:14), (11,11,11))) === nothing
    @test length(axsitr) == 27
    @test collect(axsitr) == [
        Iterators.ProductIterator((2:4, 2:4, 2:4)),
        Iterators.ProductIterator((7:9, 2:4, 2:4)),
        Iterators.ProductIterator((12:14, 2:4, 2:4)),
        Iterators.ProductIterator((2:4, 7:9, 2:4)),
        Iterators.ProductIterator((7:9, 7:9, 2:4)),
        Iterators.ProductIterator((12:14, 7:9, 2:4)),
        Iterators.ProductIterator((2:4, 12:14, 2:4)),
        Iterators.ProductIterator((7:9, 12:14, 2:4)),
        Iterators.ProductIterator((12:14, 12:14, 2:4)),
        Iterators.ProductIterator((2:4, 2:4, 7:9)),
        Iterators.ProductIterator((7:9, 2:4, 7:9)),
        Iterators.ProductIterator((12:14, 2:4, 7:9)),
        Iterators.ProductIterator((2:4, 7:9, 7:9)),
        Iterators.ProductIterator((7:9, 7:9, 7:9)),
        Iterators.ProductIterator((12:14, 7:9, 7:9)),
        Iterators.ProductIterator((2:4, 12:14, 7:9)),
        Iterators.ProductIterator((7:9, 12:14, 7:9)),
        Iterators.ProductIterator((12:14, 12:14, 7:9)),
        Iterators.ProductIterator((2:4, 2:4, 12:14)),
        Iterators.ProductIterator((7:9, 2:4, 12:14)),
        Iterators.ProductIterator((12:14, 2:4, 12:14)),
        Iterators.ProductIterator((2:4, 7:9, 12:14)),
        Iterators.ProductIterator((7:9, 7:9, 12:14)),
        Iterators.ProductIterator((12:14, 7:9, 12:14)),
        Iterators.ProductIterator((2:4, 12:14, 12:14)),
        Iterators.ProductIterator((7:9, 12:14, 12:14)),
        Iterators.ProductIterator((12:14, 12:14, 12:14))
    ]
end
