struct BVHNode{T}
    bbox::BBox{T}
    left::Int       # index of left child, -1 if leaf
    right::Int      # index of right child, -1 if leaf
    start::Int      # start index into triangle list
    count::Int      # number of triangles (only valid if leaf)
end

struct BVHModel{T,M<:GeometryBasics.Mesh}
    model::M
    nodes::Vector{BVHNode{T}}
    tri_indices::Vector{Int}
    stack::Vector{Int}
end

function buildBvh(model::M; maxLeafSize::Int=4) where {M<:GeometryBasics.Mesh}
    n = length(model)

    # Init
    tri_indices = collect(1:n)
    bmin = Vector(model[1][1])      # Just to reduce allocations
    bmax = similar(bmin)            # Just to reduce allocations

    nodes = BVHNode{eltype(model[1][1])}[]
    _build!(nodes, model, tri_indices, 1, n, maxLeafSize, bmin, bmax)
    return BVHModel(model, nodes, tri_indices, Int[])
end

function _build!(nodes, model, tri_indices, start, stop, maxLeafSize, bmin, bmax)
    count = stop - start + 1

    # compute bbox for [start:stop]
    @inbounds for i in 1:3
        bmin[i] = model[tri_indices[start]][1][i]
        bmax[i] = bmin[i]
    end
    for i in start:stop
        tri = model[tri_indices[i]]
        for j in 1:3
            bmin[j] = min(bmin[j], tri[1][j], tri[2][j], tri[3][j])
            bmax[j] = max(bmax[j], tri[1][j], tri[2][j], tri[3][j])
        end
    end
    node_bbox = BBox(bmin, bmax)
    node_idx = length(nodes) + 1

    if count <= maxLeafSize
        push!(nodes, BVHNode(node_bbox, -1, -1, start, count))
        return node_idx
    end

    # choose split axis (longest axis)
    diagx = node_bbox.maxx - node_bbox.minx
    diagy = node_bbox.maxy - node_bbox.miny
    diagz = node_bbox.maxz - node_bbox.minz
    axis = findmax((diagx, diagy, diagz))[2]

    # partition tri_indices[start:stop] by centroid along axis
    mid = (start + stop) >>> 1
    v = view(tri_indices, start:stop)
    sort!(v; by=i -> begin
        v1, v2, v3 = model[i]
        (v1[axis] + v2[axis] + v3[axis]) # * onethird (multiplication by 1/3 not necessary for sorting)
    end)

    # push placeholder parent (pre-order)
    push!(nodes, BVHNode(node_bbox, -1, -1, -1, 0))

    # build children
    left = _build!(nodes, model, tri_indices, start, mid, maxLeafSize, bmin, bmax)
    right = _build!(nodes, model, tri_indices, mid+1, stop, maxLeafSize, bmin, bmax)

    # overwrite parent with correct child indices
    nodes[node_idx] = BVHNode(node_bbox, left, right, -1, 0)
    return node_idx
end

function intersect!(ray::Ray{T}, bvh::BVHModel{T,M}, anyHit=false) where {T,M<:GeometryBasics.Mesh}
    resetRay!(ray)
    invdirRay!(ray)
    stack = bvh.stack
    empty!(stack)
    push!(stack, 1)     # start at root
    tInf = T(Inf)
    while !isempty(stack)
        node_idx = pop!(stack)
        node = bvh.nodes[node_idx]

        if intersect(ray, node.bbox)
            if node.left < 0   # leaf
                @inbounds for i in node.start:(node.start + node.count - 1)
                    idx = bvh.tri_indices[i]
                    if idx != ray.idxSkip
                        tri = bvh.model[idx]
                        intersect!(ray, tri, idx)
                        if anyHit && ray.t < tInf
                            break
                        end
                    end
                end
            else               # internal
                push!(stack, node.left)
                push!(stack, node.right)
            end
        end
    end
    return ray.t < tInf
end
