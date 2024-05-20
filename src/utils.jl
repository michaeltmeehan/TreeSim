"""
    CantorPair(x, y)


    CantorPair(x,y)
    Calculate the Cantor Pair mapping of `x` and `y`.
"""
function CantorPair(x, y)
    y == 0 && return 0
    return trunc(Int, (x^2 + 3 * x + 2 * x * y + y + y^2) / 2)
end