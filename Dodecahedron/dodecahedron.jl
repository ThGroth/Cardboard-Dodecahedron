#
#
#@author Laurent Bartholdi
#
#
const recompute = false # true to compute data, false to read it hard-coded
using Colors
# group preserving the form w^2-x^2-y^2-z^2
typealias O31 Array{Float64,2}

const O31form = [1.0 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]::O31

hyperbolicdist(x,y) = acosh(x⋅(O31form*y))
hyperbolicnorm(x) = acosh(x[1]) # really hyperbolicdist(x,[1,0,0,0])
kleinprojection(x) = x[2:4]/x[1]

# constants
const size_alt5 = 60
const num_gens = 12
const num_vertices = 20
const num_edges = 30

if recompute

# 1.618...
τ = (1+√5)/2
ɛ = eps(1000.0)

# given a quaternion w+ix+jy+kz
# compute the corresponding rotation in so31 fixing [1,0,0,0]
function quaternion_to_so31(w,x,y,z)
    [w^2+x^2+y^2+z^2 0 0 0;
     0 w^2+x^2-y^2-z^2 2*x*z-2*w*y -2*w*z-2*x*y;
     0 2*w*y+2*x*z w^2-x^2-y^2+z^2 2*x*w-2*y*z;
     0 2*w*z-2*x*y -2*x*w-2*y*z w^2-x^2+y^2-z^2]::O31/(w^2+x^2+y^2+z^2)
end

# the transformation shifting in direction v.
# computed as exp([0 v; v' 0]) in block matrix form
function boost(v::Array{Float64,1})
    λ = norm(v)
    w = v/λ
    m = zeros(4,4)
    m[1,1] = cosh(λ)
    m[2:4,1] = w*sinh(λ)
    m[1,2:4] = w*sinh(λ)
    m[2:4,2:4] = eye(3)+(cosh(λ)-1)*w*w'
    m
end

# the elementary reflection fixing [1,0,0,0] in a plane perpendicular to v
function reflection(v::Array{Float64,1})
    m = zeros(4,4)
    m[1,1] = 1
    m[2:4,2] = v/norm(v)
    m[2:4,3:4] = nullspace(v')
    m*diagm([1,-1,1,1])*m'
end

# the sign of a permutation
function signature(perm::Array{Int,1})
    m = zeros(Int,length(perm),length(perm))
    for i=1:length(perm)
        m[i,perm[i]] = 1
    end
    sign(det(m))
end

# order of a group element, in {1,2,3,5}
function order(mat::O31)
    for i=1:5
        norm(mat^i-alt5[1])<ɛ && return i
    end
    0 # infinity
end

# the alternating group, with 60 elements
alt5 = O31[]
push!(alt5,quaternion_to_so31(1.0,0,0,0))
push!(alt5,quaternion_to_so31(0,1.0,0,0))
push!(alt5,quaternion_to_so31(0,0,1.0,0))
push!(alt5,quaternion_to_so31(0,0,0,1.0))
for i=-1:2:1, j=-1:2:1, k=-1:2:1
    push!(alt5,quaternion_to_so31(1/2,i/2,j/2,k/2))
end
for i=-1:2:1, j=-1:2:1
    v=[i/τ,1,j*τ,0]/2
    for perm=permutations([1,2,3,4])
        signature(perm)==1 && push!(alt5,quaternion_to_so31(v[perm]...))
    end
end

# the 12 coxeter generators
generators = Array(O31,num_gens)
generators[1] = [τ+1 0 -sqrt(2τ+1) -sqrt(τ);
                 0 1 0 0;
                 sqrt(2τ+1) 0 -τ -τ;
                 sqrt(τ) 0 -τ 0]
numgens = 1 # lazy, just compute orbit under alt5
for i=2:size_alt5
    g = alt5[i]^-1*generators[1]*alt5[i]
    keep = true
    for j=1:numgens
        if norm(g-generators[j])<ɛ keep=false; break end
    end
    keep && (numgens += 1; generators[numgens] = g)
end

# the 20 vertices of the dodecahedron
vertices = Array(Array{Float64,1},num_vertices)
for i=1:num_gens
    v = [τ^2/√2,0,0,0]
    v[2+((i+1)%3)] = (((i-1)%6)≥3?-1:1)*sqrt(1/2τ)
    v[2+((i+2)%3)] = (i>6?-1:1)*sqrt(τ+1/2)
    vertices[i] = v
end
numvert = 13
for i=-1:2:1, j=-1:2:1, k=-1:2:1
    v = [τ^2/√2,i*sqrt(τ/2),j*sqrt(τ/2),k*sqrt(τ/2)]
    vertices[numvert] = v
    numvert += 1
end

# the 30 edges of the dodecahedron
edges = Array(Array{Float64,2},num_edges)
edges[1] = [vertices[1] vertices[4]]
numedge = 1 # lazy, just compute orbit under alt5
for i=2:size_alt5
    e = alt5[i]*edges[1]
    keep = true
    for j=1:numedge
        if norm(e-edges[j])<ɛ || norm(e[:,[2,1]]-edges[j])<ɛ
            keep=false; break
        end
    end
    keep && (numedge += 1; edges[numedge] = e)
end

# the coxeter matrix. Only store if generators commute
preccommute = Bool[i≤j && norm(generators[i]*generators[j]-generators[j]*generators[i])<ɛ for i=1:num_gens, j=1:num_gens]
commute = preccommute|preccommute'

# which coxeter generator fixes which edge/vertex, in format
# fixesvertex(generators[i],vertices[j])
# fixesedge(generators[i],vertices[j])
fixesvertex = Bool[norm(generators[g]*vertices[j]-vertices[j])<ɛ for g=1:num_gens, j=1:num_vertices]
fixesedge = Bool[norm(generators[g]*edges[j]-edges[j])<ɛ for g=1:num_gens, j=1:num_edges]

else # hard-coded data
preccommute = Bool[true false false true false true true false false true true false
     false true true false true false false true true false false true
     false false true false false false false false true true true true
     false false false true true true true true false false false false
     false false false false true true false true true false false false
     false false false false false true false false true false true false
     false false false false false false true true false true false true
     false false false false false false false true false false false true
     false false false false false false false false true false true false
     false false false false false false false false false true true true
     false false false false false false false false false false true false
     false false false false false false false false false false false true]
commute = Bool[true false false true false true true false false true true false
     false true true false true false false true true false false true
     false true true false false false false false true true true true
     true false false true true true true true false false false false
     false true false true true true false true true false false false
     true false false true true true false false true false true false
     true false false true false false true true false true false true
     false true false true true false true true false false false true
     false true true false true true false false true false true false
     true false true false false false true false false true true true
     true false true false false true false false true true true false
     false true true false false false true true false true false true]
fixesvertex = Bool[false true true false true false false false false false false false false false false true false false false true
     false false false false false false false true false false true true true false false false true false false false
     false false false false false true false true false false true false false true false false false true false false
     false true false false true false false false true false false false false false true false false false true false
     false false false false false false false false true true false true true false true false false false false false
     false false false false true false true false false true false false false false true true false false false false
     true true false true false false false false false false false false false false false false false false true true
     false false false true false false false false true false false true false false false false true false true false
     false false false false false false true false false true true false true true false false false false false false
     true false true false false true false false false false false false false false false false false true false true
     false false true false false true true false false false false false false true false true false false false false
     true false false true false false false true false false false false false false false false true true false false]
fixesedge = Bool[false false false false true false true false false false false false false false false false false true false false false false true false true false false false false false
     false false false true false false false true false false false true false false false false false false true false false false false false false true false false false false
     false false false true false false false false false false false false true false false false true false false false true false false true false false false false false false
     false false false false true false false false false true false false false false true false false false false true false true false false false false false false false false
     false false true false false false false false false true false false false false false true false false true false false false false false false false true false false false
     false true false false false false false false false false false false false false false false false false false false false true false false true false true false false true
     true false false false false false false false false false false false false false false false false false false true false false true false false false false true true false
     false false true false false false false false false false true false false false true false false false false false false false false false false true false false true false
     false true false false false false false false true false false true false false false true true false false false false false false false false false false false false false
     false false false false false true false false false false false false false true false false false true false false false false false true false false false true false false
     false false false false false true true false true false false false false false false false false false false false true false false false false false false false false true
     true false false false false false false true false false true false true true false false false false false false false false false false false false false false false false]

generators = O31[
[2.618033988749895 0.0 -2.0581710272714924 -1.272019649514069
 0.0 1.0 0.0 0.0
 2.0581710272714924 0.0 -1.618033988749895 -1.618033988749895
 1.272019649514069 0.0 -1.618033988749895 0.0],

[2.618033988749895 0.0 2.0581710272714924 1.272019649514069
 0.0 1.0 0.0 0.0
 -2.0581710272714924 0.0 -1.618033988749895 -1.618033988749895
 -1.272019649514069 0.0 -1.618033988749895 0.0],

[2.618033988749895 0.0 2.0581710272714924 -1.272019649514069
 0.0 1.0 0.0 0.0
 -2.0581710272714924 0.0 -1.618033988749895 1.618033988749895
 1.272019649514069 0.0 1.618033988749895 0.0],

[2.618033988749895 0.0 -2.0581710272714924 1.272019649514069
 0.0 1.0 0.0 0.0
 2.0581710272714924 0.0 -1.618033988749895 1.618033988749895
 -1.272019649514069 0.0 1.618033988749895 0.0],

[2.618033988749895 1.272019649514069 0.0 2.0581710272714924
 -1.272019649514069 0.0 0.0 -1.618033988749895
 0.0 0.0 1.0 0.0
 -2.0581710272714924 -1.618033988749895 0.0 -1.618033988749895],

[2.618033988749895 2.0581710272714924 -1.272019649514069 0.0
 -2.0581710272714924 -1.618033988749895 1.618033988749895 0.0
 1.272019649514069 1.618033988749895 0.0 0.0
 0.0 0.0 0.0 1.0],

[2.618033988749895 -2.0581710272714924 -1.272019649514069 0.0
 2.0581710272714924 -1.618033988749895 -1.618033988749895 0.0
 1.272019649514069 -1.618033988749895 0.0 0.0
 0.0 0.0 0.0 1.0],

[2.618033988749895 -1.272019649514069 0.0 2.0581710272714924
 1.272019649514069 0.0 0.0 1.618033988749895
 0.0 0.0 1.0 0.0
 -2.0581710272714924 1.618033988749895 0.0 -1.618033988749895],

[2.618033988749895 2.0581710272714924 1.272019649514069 0.0
 -2.0581710272714924 -1.618033988749895 -1.618033988749895 0.0
 -1.272019649514069 -1.618033988749895 0.0 0.0
 0.0 0.0 0.0 1.0],

[2.618033988749895 -1.272019649514069 0.0 -2.0581710272714924
 1.272019649514069 0.0 0.0 -1.618033988749895
 0.0 0.0 1.0 0.0
 2.0581710272714924 -1.618033988749895 0.0 -1.618033988749895],

[2.618033988749895 1.272019649514069 0.0 -2.0581710272714924
 -1.272019649514069 0.0 0.0 1.618033988749895
 0.0 0.0 1.0 0.0
 2.0581710272714924 1.618033988749895 0.0 -1.618033988749895],

[2.618033988749895 -2.0581710272714924 1.272019649514069 0.0
 2.0581710272714924 -1.618033988749895 1.618033988749895 0.0
 -1.272019649514069 1.618033988749895 0.0 0.0
 0.0 0.0 0.0 1.0]]
vertices = Array{Float64,1}[[1.851229586821916,1.455346690225355,0.0,0.5558929702514211],
  [1.851229586821916,0.5558929702514211,1.455346690225355,0.0],
  [1.851229586821916,0.0,0.5558929702514211,1.455346690225355],
  [1.851229586821916,1.455346690225355,0.0,-0.5558929702514211],
  [1.851229586821916,-0.5558929702514211,1.455346690225355,0.0],
  [1.851229586821916,0.0,-0.5558929702514211,1.455346690225355],
  [1.851229586821916,-1.455346690225355,0.0,0.5558929702514211],
  [1.851229586821916,0.5558929702514211,-1.455346690225355,0.0],
  [1.851229586821916,0.0,0.5558929702514211,-1.455346690225355],
  [1.851229586821916,-1.455346690225355,0.0,-0.5558929702514211],
  [1.851229586821916,-0.5558929702514211,-1.455346690225355,0.0],
  [1.851229586821916,0.0,-0.5558929702514211,-1.455346690225355],
  [1.851229586821916,-0.8994537199739336,-0.8994537199739336,-0.8994537199739336],
  [1.851229586821916,-0.8994537199739336,-0.8994537199739336,0.8994537199739336],
  [1.851229586821916,-0.8994537199739336,0.8994537199739336,-0.8994537199739336],
  [1.851229586821916,-0.8994537199739336,0.8994537199739336,0.8994537199739336],
  [1.851229586821916,0.8994537199739336,-0.8994537199739336,-0.8994537199739336],
  [1.851229586821916,0.8994537199739336,-0.8994537199739336,0.8994537199739336],
  [1.851229586821916,0.8994537199739336,0.8994537199739336,-0.8994537199739336],
  [1.851229586821916,0.8994537199739336,0.8994537199739336,0.8994537199739336]]
edges = Array{Float64,2}[
[1.851229586821916 1.851229586821916
 1.455346690225355 1.455346690225355
 0.0 0.0
 0.5558929702514211 -0.5558929702514211],

[1.851229586821916 1.851229586821916
 -1.455346690225355 -1.455346690225355
 0.0 0.0
 0.5558929702514211 -0.5558929702514211],

[1.851229586821916 1.851229586821916
 0.0 0.0
 -0.5558929702514211 0.5558929702514211
 -1.455346690225355 -1.455346690225355],

[1.851229586821916 1.851229586821916
 -0.5558929702514211 0.5558929702514211
 -1.455346690225355 -1.455346690225355
 0.0 0.0],

[1.851229586821916 1.851229586821916
 0.5558929702514211 -0.5558929702514211
 1.455346690225355 1.455346690225355
 0.0 0.0],

[1.851229586821916 1.851229586821916
 0.0 0.0
 -0.5558929702514211 0.5558929702514211
 1.455346690225355 1.455346690225355],

[1.851229586821916 1.851229586821916
 -6.221894620906635e-17 -0.8994537199739336
 0.5558929702514211 0.8994537199739335
 1.455346690225355 0.8994537199739339],

[1.851229586821916 1.851229586821916
 0.8994537199739336 0.5558929702514213
 -0.8994537199739339 -1.455346690225355
 -0.8994537199739336 -6.7077949778085215e-18],

[1.851229586821916 1.851229586821916
 -1.455346690225355 -0.8994537199739339
 -6.7077949778085215e-18 -0.8994537199739336
 0.5558929702514211 0.8994537199739335],

[1.851229586821916 1.851229586821916
 -6.221894620906635e-17 -0.8994537199739336
 0.5558929702514211 0.8994537199739335
 -1.455346690225355 -0.8994537199739339],

[1.851229586821916 1.851229586821916
 0.8994537199739339 1.455346690225355
 -0.8994537199739336 -6.7077949778085215e-18
 -0.8994537199739335 -0.5558929702514211],

[1.851229586821916 1.851229586821916
 -0.5558929702514213 -0.8994537199739336
 -1.455346690225355 -0.8994537199739339
 -6.7077949778085215e-18 -0.8994537199739336],

[1.851229586821916 1.851229586821916
 0.8994537199739336 0.5558929702514213
 -0.8994537199739339 -1.455346690225355
 0.8994537199739336 6.7077949778085215e-18],

[1.851229586821916 1.851229586821916
 0.8994537199739339 1.455346690225355
 -0.8994537199739336 -6.7077949778085215e-18
 0.8994537199739335 0.5558929702514211],

[1.851229586821916 1.851229586821916
 0.8994537199739336 6.221894620906635e-17
 0.8994537199739335 0.5558929702514211
 -0.8994537199739339 -1.455346690225355],

[1.851229586821916 1.851229586821916
 -1.455346690225355 -0.8994537199739339
 -6.7077949778085215e-18 -0.8994537199739336
 -0.5558929702514211 -0.8994537199739335],

[1.851229586821916 1.851229586821916
 -0.5558929702514213 -0.8994537199739336
 -1.455346690225355 -0.8994537199739339
 6.7077949778085215e-18 0.8994537199739336],

[1.851229586821916 1.851229586821916
 0.8994537199739336 6.221894620906635e-17
 0.8994537199739335 0.5558929702514211
 0.8994537199739339 1.455346690225355],

[1.851229586821916 1.851229586821916
 -0.8994537199739336 -6.221894620906635e-17
 -0.8994537199739335 -0.5558929702514211
 -0.8994537199739339 -1.455346690225355],

[1.851229586821916 1.851229586821916
 0.8994537199739336 0.5558929702514213
 0.8994537199739339 1.455346690225355
 -0.8994537199739336 -6.7077949778085215e-18],

[1.851229586821916 1.851229586821916
 -0.8994537199739336 -6.221894620906635e-17
 -0.8994537199739335 -0.5558929702514211
 0.8994537199739339 1.455346690225355],

[1.851229586821916 1.851229586821916
 -0.5558929702514213 -0.8994537199739336
 1.455346690225355 0.8994537199739339
 -6.7077949778085215e-18 -0.8994537199739336],

[1.851229586821916 1.851229586821916
 0.8994537199739336 0.5558929702514213
 0.8994537199739339 1.455346690225355
 0.8994537199739336 6.7077949778085215e-18],

[1.851229586821916 1.851229586821916
 6.221894620906635e-17 0.8994537199739336
 -0.5558929702514211 -0.8994537199739335
 1.455346690225355 0.8994537199739339],

[1.851229586821916 1.851229586821916
 -0.5558929702514213 -0.8994537199739336
 1.455346690225355 0.8994537199739339
 6.7077949778085215e-18 0.8994537199739336],

[1.851229586821916 1.851229586821916
 6.221894620906635e-17 0.8994537199739336
 -0.5558929702514211 -0.8994537199739335
 -1.455346690225355 -0.8994537199739339],

[1.851229586821916 1.851229586821916
 -0.8994537199739339 -1.455346690225355
 0.8994537199739336 6.7077949778085215e-18
 -0.8994537199739335 -0.5558929702514211],

[1.851229586821916 1.851229586821916
 1.455346690225355 0.8994537199739339
 6.7077949778085215e-18 0.8994537199739336
 0.5558929702514211 0.8994537199739335],

[1.851229586821916 1.851229586821916
 1.455346690225355 0.8994537199739339
 6.7077949778085215e-18 0.8994537199739336
 -0.5558929702514211 -0.8994537199739335],

[1.851229586821916 1.851229586821916
 -0.8994537199739339 -1.455346690225355
 0.8994537199739336 6.7077949778085215e-18
 0.8994537199739335 0.5558929702514211]]
end
#=
function vertexexplore(maxnorm::Float64,callback::Function)
    function recur(lastgens::BitVector,init::Int,point::Array{Float64,1})
        # lastgens are the generators we last saw and that potentially
        #   could cancel with a new generator
        # init is the vertex we started with, or 0 if we left it
        # point is a point on the hyperboloid
        
        hyperbolicnorm(point)<maxnorm || return
        callback(point)
        for g=1:num_gens
            lastgens[g] && continue
            init>0 && fixesvertex[g,init] && continue
            recur((lastgens&commute[:,g])|preccommute[:,g],0,generators[g]*point)
        end
    end
    for i=1:num_vertices
        recur(BitVector(num_gens),i,vertices[i])
    end
end
=#
function vertexexplore(maxnorm::Float64,callback::Function)
    function recur(lastgens::BitVector,init::Int,point::Array{Float64,1})
        # lastgens are the generators we last saw and that potentially
        #   could cancel with a new generator
        # init is the vertex we started with, or 0 if we left it
        # point is a point on the hyperboloid
        hyperbolicnorm(point)<maxnorm || return
        callback(point)
        for g=1:num_gens
            lastgens[g] && continue
            init>0 && fixesvertex[g,init] && continue
            recur((lastgens&commute[:,g])|preccommute[:,g],0,generators[g]*point)
        end
    end
    for i=1:num_vertices
        recur(BitVector(num_gens),i,vertices[i])
    end
end
#=
function edgeexplore(maxnorm::Float64,callback::Function)
    function recur(lastgens::BitVector,init::Int,edge::Array{Float64,2})
        # lastgens are the generators we last saw and that potentially
        #   could cancel with a new generator
        # init is the edge we started with, or 0 if we left it
        # point is a point on the hyperboloid
        (hyperbolicnorm(edge[:,1])<maxnorm && hyperbolicnorm(edge[:,2])<maxnorm) || return
        callback(edge[:,1],edge[:,2])
        for g=1:num_gens
            lastgens[g] && continue
            init>0 && fixesedge[g,init] && continue
            recur((lastgens&commute[:,g])|preccommute[:,g],0,generators[g]*edge)
        end
    end
    for i=1:num_edges
        recur(BitVector(num_gens),i,edges[i])
    end
end
=#
#Version for coloring all edges corresponding to the same generator with the same color
# function edgeexplore(maxnorm::Float64,callback::Function)
#     function recur(lastgens::BitVector,init::Int,edge::Array{Float64,2},gennum::Int)
#         # lastgens are the generators we last saw and that potentially
#         #   could cancel with a new generator
#         # init is the edge we started with, or 0 if we left it
#         # point is a point on the hyperboloid
#         (hyperbolicnorm(edge[:,1])<maxnorm && hyperbolicnorm(edge[:,2])<maxnorm) || return
#         callback(edge[:,1],edge[:,2],gennum)
#         for g=1:num_gens
#             lastgens[g] && continue
#             init>0 && fixesedge[g,init] && continue
#             recur((lastgens&commute[:,g])|preccommute[:,g],0,generators[g]*edge,gennum)
#         end
#     end
#     for i=1:num_edges
#         recur(BitVector(num_gens),i,edges[i],i)
#     end
# end
function edgeexplore(maxnorm::Float64,callback::Function)
    function recur(lastgens::BitVector,init::Int,edge::Array{Float64,2})
        # lastgens are the generators we last saw and that potentially
        #   could cancel with a new generator
        # init is the edge we started with, or 0 if we left it
        # point is a point on the hyperboloid
        (hyperbolicnorm(edge[:,1])<maxnorm && hyperbolicnorm(edge[:,2])<maxnorm) || return
        callback(edge[:,1],edge[:,2])
        for g=1:num_gens
            lastgens[g] && continue
            init>0 && fixesedge[g,init] && continue
            recur((lastgens&commute[:,g])|preccommute[:,g],0,generators[g]*edge)
        end
    end
    for i=1:num_edges
        recur(BitVector(num_gens),i,edges[i])
    end
end

# vertexcallback takes a vector on the hyperboloid as argument
# edgecallback takes two vectors (start and end of edge) as argument
function explore(maxnorm::Float64,vertexcallback::Function,edgecallback::Function)
    vertexexplore(maxnorm,vertexcallback)
    edgeexplore(maxnorm,edgecallback)
end

# sample use

#=
#f = open("dode40vec4","w")
#explore(4.0,
#        v->println(f,join(v," ")),
#        (v,w)->println(f,join(v," ")," ",join(w," ")))
#close(f)
=#

# f = open("dode2.txt","w")
# explore(2.0,
#         v->println(f,join(kleinprojection(v)," ")),
#         (v,w)->println(f,join(kleinprojection(v)," ")," ",join(kleinprojection(w)," ")))
# close(f)
function col(c)
    return hex(round(Int,c.r*256))*hex(round(Int,c.g*256))*hex(round(Int,c.b*256))
end 
bound = 4.0
f = open("dode32.js","w")

println(f,"DodeLines  = [];\n\nvar material = new THREE.LineBasicMaterial({color: 0x0000ff});\n")

edgeexplore(bound,
    (v,w)->println(f,"var material = new THREE.LineBasicMaterial({color: 0x",
        hex(HSL(9.0,
                1.0-(hyperbolicnorm(v)/(1.2*bound))^2,
                0.5*(1.0-(hyperbolicnorm(v)/bound)^2))),"});\n",
        "var geometry = new THREE.Geometry();\ngeometry.vertices.push(\n\tnew THREE.Vector3(",
        join(10*kleinprojection(v),","),"),\n\tnew THREE.Vector3(",
        join(10*kleinprojection(w),","),")\n);\n",
        "var line = new THREE.Line( geometry, material );\n",
        "DodeLines.push( line );\n"))
close(f)


# edgeexplore(bound,
#     (v,w)->println(f,join(kleinprojection(v)," "),
#         " ",join(kleinprojection(w)," "),
#         " ",join(col(convert(RGB, HSL(9.0,
#                                     1.0-(hyperbolicnorm(v)/(1.2*bound))^2,
#                                     0.5*(1.0-(hyperbolicnorm(v)/bound)^2))))," ")))
# close(f)

# the transformation shifting in direction v.
# computed as exp([0 v; v' 0]) in block matrix form
function boost(v::Array{Float64,1})
    λ = norm(v)
    w = v/λ
    m = zeros(4,4)
    m[1,1] = cosh(λ)
    m[2:4,1] = w*sinh(λ)
    m[1,2:4] = w*sinh(λ)
    m[2:4,2:4] = eye(3)+(cosh(λ)-1)*w*w'
    m
end

#shift = [0.0,0.0,3.0]
#m=boost(shift)

# f = open("dode2.txt","w")
# explore(4.0,
#         v->print(f,"point3d(",kleinprojection(v),")+"),
#         (v,w)->print(f,"line3d([",kleinprojection(v),",",kleinprojection(w),"],color=hue(",i+30,"/60))+"))
# close(f)



#f = open("dodecahedron.data","w")
#explore(3.0,
#        v->println(f,kleinprojection(v)),
#        (v,w)->println(f,kleinprojection(v),"-->",kleinprojection(w)))
#close(f)

# to move a set of points in a direction (given by a vector v),
# multiply them on the left by boost(v)

# sample points
