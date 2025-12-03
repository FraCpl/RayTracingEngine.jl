# a = b × c
@inline function cross!(a, b, c)
    a[1] = -b[3]*c[2] + b[2]*c[3]
    a[2] = +b[3]*c[1] - b[1]*c[3]
    a[3] = -b[2]*c[1] + b[1]*c[2]
    return nothing
end

@inline function dot3(a, b)
    ax, ay, az = a
    bx, by, bz = b
    return ax*bx + ay*by + az*bz
end
