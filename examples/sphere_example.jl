using DDA

using plot

# incident field
E0 = 1
k = [0, 0, 1]
e = [1, 0, 0]
E_inc = E_inc()



# grid and target
s = Sprhere(30)
grid = CubicGrid(-20:1:20,-20:1:20,-2:1:2,)


dipoles = get_dipoles(g, s)

# material model

# solve
lambda = 550.
# lambda = 400:10:600
solve(dipoles, E_inc)

??material