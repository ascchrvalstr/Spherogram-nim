import std/random
import std/options
import std/strformat

{.compile: "planarmap/PMconjugation.c".}
{.compile: "planarmap/PMdef.c".}
{.compile: "planarmap/PMenlight.c".}
{.compile: "planarmap/PMextract.c".}
{.compile: "planarmap/PMplanmap.c".}
{.compile: "planarmap/stats.c".}

type
    pmSize* = object
        m*, b*: cchar     # map and basic map type
        e*, v*, f*: clong # edges, vertices, faces
        r*, g*, d*: clong # red and black vertices, max degree
        t*: clong         # tolerance on e
        dgArr*: ptr clong # pt on vertex list
    pmMethod* = object
        core*, pic*: cchar
        seed*: clong
        verbose*: cchar
    pmMemory* = object
        dTree*: cchar
        sTree*, rTree*, gTree*, sWrd*, sEdge*, sVtx*, sLeaf*: clong
    pm_vertex* = object
        root*: ptr pm_edge
        next*: ptr pm_vertex
        mark*: clong
        type0*: cshort
        label*: clong
    pm_edge* = object
        from0*: ptr pm_vertex
        face*: ptr pm_vertex
        prev*: ptr pm_edge
        next*: ptr pm_edge
        oppo*: ptr pm_edge
        mark*: clong
        type0*: cshort
        label*: clong
    pmMap* = object
        root*: ptr pm_edge
        e*, v*, f*, i*: clong

proc set_pmRandom_callback*(function: (proc(n: clong): clong)): void {.importc.}

proc pmMemoryInit*(S: ptr pmSize, Meth: ptr pmMethod, M: ptr pmMemory): cint {.importc.}
proc pmPlanMap*(S: ptr pmSize, Meth: ptr pmMethod, M: ptr pmMemory, Map: ptr pmMap): cint {.importc.}
proc pmFreeMap*(Map: ptr pmMap): cint {.importc.}

proc pmStatGauss*(Map: ptr pmMap): clong {.importc.}

proc randrange_callback*(n: clong): clong =
    return rand(1 .. n)

proc raw_random_map*(num_vertices: int, edge_connectivity: int = 4, num_link_comps: int = 0, max_tries: int = 100): Option[seq[(int, seq[int])]] =
    discard """
    Use Gilles Schaeffer's "Planarmap program" to generate
    a random 4-valent planar graph with the given number
    of vertices.

    The "edge_connectivity" parameter can be any of 2, 4, or 6.
    Recall that a graph is k-edge-connected if removing any set
    of *less than* k edges disconnects the graph. In particular,
    for 4-valent graphs, being 2-connected is just the same as
    being connected. In particular, a 2-connected graph can (and
    frequently do) have looped edges.

    Under the hood, it uses Nim's pseudo-random number
    generator.
    """
    set_pmRandom_callback(randrange_callback)

    var size: pmSize
    var method0: pmMethod
    var memory: pmMemory
    var the_map: pmMap
    var edge: ptr pm_edge
    var vert: ptr pm_vertex

    if edge_connectivity == 2:
        (size.m, size.b) = (char(4), char(4))
    elif edge_connectivity == 4:
        (size.m, size.b) = (char(5), char(5))
    elif edge_connectivity == 6:
        (size.m, size.b) = (char(6), char(5))
    else:
        raise newException(ValueError, &"Invalid edge_connectivity parameter {edge_connectivity}")

    size.v = num_vertices
    (size.e, size.f, size.r, size.g, size.d) = (0, 0, 0, 0, 0)
    size.t = -1
    size.dgArr = nil

    (method0.core, method0.pic) = (char(0), char(0))
    method0.verbose = char(0)
    discard pmMemoryInit(addr size, addr method0, addr memory)
    discard pmPlanMap(addr size, addr method0, addr memory, addr the_map)
    if num_link_comps > 0:
        var tries = 0
        while pmStatGauss(addr the_map) != num_link_comps and tries < max_tries:
            discard pmFreeMap(addr the_map)
            discard pmPlanMap(addr size, addr method0, addr memory, addr the_map)
            tries += 1

        if tries == max_tries:
            return none(seq[(int, seq[int])])

    var ans = newSeq[(int, seq[int])]()
    vert = the_map.root.from0
    while vert != nil:
        var edges_at_vert = newSeq[int]()
        edge = vert.root
        while edge != vert.root.prev:
            edges_at_vert.add(int(edge.label))
            edge = edge.next
        edges_at_vert.add(edge.label)
        ans.add((int(vert.label), edges_at_vert))
        vert = vert.next
    discard pmFreeMap(addr the_map)
    return some(ans)
