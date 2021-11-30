using DDA, Test

# # 3D Cartesian grid examples
#
# grid = CartesianGrid{Float64,3}([0, 0, 0], [.1,.1,.1], (10, 10, 10))
# grid = CartesianGrid([0, 0, 0], [.1,.1,.1], (10, 10, 10))
# grid = CartesianGrid([0, 0, 0], [.1,.1,.1], [10, 10, 10])
# grid = CartesianGrid([0., 0., 0.], [.1,.1,.1], (10, 10, 10))
#
# grid = CartesianGrid([0, 0, 0], [1, 1, 1], (10, 10, 10))
# grid = CartesianGrid(10, 10, 10)


# using StaticArrays
# g = CartesianGrid(SVector(0., 0, 0), SVector(1., 1, 1), (10, 10, 10))
#
# grid = CartesianGrid([0, 0, 0], [.1,.1,.1], (10, 10, 10))
# eltype(grid)
#
# grid = CartesianGrid([0, 0, 0], [.1,.1,.1], (10, 10, 10))
# grid = CartesianGrid([0, 0, 0], [.1,.1,.1], [10, 10, 10])
# grid = CartesianGrid([0., 0., 0.], [.1,.1,.1], (10, 10, 10))
#
# grid = CartesianGrid([0, 0, 0], [1, 1, 1], (10, 10, 10))
# grid = CartesianGrid(10, 10, 10)
#
#
# g = CartesianGrid(100, 100, 50)
# DDA.center(g)
#
#
# g = CartesianGrid([10.0, 20.0], [1.0, 1.0], (100, 1000))
# g = CartesianGrid([10.0, 20.0], [1.0, 1.0], [100, 1000])
# # g = CartesianGrid((10.,20.), (1,1), (100,100))
#
# DDA.center(g)
#
# g = CartesianGrid([10.0, 20.0, 30], [1.0, 1.0, 2], [101, 101, 51])
# DDA.center(g)
#
# maximum(g)
# extrema(g)
# DDA.width(g)
#
# g[CartesianIndex(1,1,1)]
# g[CartesianIndex(1,1,1,1,1)]
# g[CartesianIndex(1,1,100000)]
# g[1,1,1]
# g[1,1,5000]
#
# view(g, 1,:,:)
#
#
# g[1,:,1]
# collect(g[1,:,1])
# @inbounds  g[1,1:300000,1]
# @inbounds  collect(g[1,1:300000,1])
#
# c = 0
# for ind in CartesianIndices(g)
#     # @show ind g[ind]
#     # break
#     c+=1
# end
# c
# for coord in g
#     @show coord typeof(coord)
#     break
# end
#
#
# eltype(g)
#
# CartesianGrid{Float64}(10, 10, 10)
# CartesianGrid(10, 10, 10)
#
#
# CartesianGrid(100, 100, 50)
# CartesianGrid([10.,20.], [1.,1.], (100, 100))