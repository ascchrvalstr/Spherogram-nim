discard """
Generating random knots starting with randomly generated planar maps.

Code written as part of a 2013 Illinois Geometry Lab project by

Nathan Dunfield, Alexander Ehrenberg, Shiladitya Bhattacharyya, and Dongming Lei

Details to hopefully appear in some paper or other. I wouldn't hold
your breath, though.
"""
import std/options
import std/strformat
import std/random
import std/sequtils
import std/enumerate

import planarmap
import links

type LinkGenerationError* = object of CatchableError

proc random_link_internal*[T](num_crossings: int, edge_conn_param: int = 4, num_link_comps: int = 0, max_tries: int = 100): Link[T] =
    discard """
    Returns a dictionary of endpoints of edges in the form:

    (signed edge) -> (vertex, position)
    """
    if num_crossings < 0:
        raise newException(ValueError, &"the number of crossings {num_crossings} cannot be negative")
    if num_link_comps < 0:
        raise newException(ValueError, &"the number of link components {num_link_comps} cannot be negative")
    if max_tries <= 0:
        raise newException(ValueError, &"the maximum number of attempts {max_tries} cannot be negative or zero")
    if num_crossings == 0:
        if num_link_comps > 0:
            raise newException(ValueError, &"a zero-crossing link diagram is empty, so it cannot admit {num_link_comps} link components")
        return link_from_PD_code_with_extra_info[T](@[])
    if num_crossings == 1:
        if num_link_comps > 1:
            raise newException(ValueError, &"a one-crossing link diagram must be an unknot after a reverse Reidemeister I move, so it must have one link component, so it cannot have {num_link_comps} link components")
        # there are exactly two possibilities depending on whether the opposite of strand 0 is strand 1 or strand 3
        if rand(0..1) == 0:
            return link_from_PD_code_with_extra_info[T](@[[0, 0, 1, 1]])
        else:
            return link_from_PD_code_with_extra_info[T](@[[0, 1, 1, 0]])
    var data0 = raw_random_map(num_crossings, edge_conn_param, num_link_comps, max_tries)
    if data0.isNone:
        raise newException(LinkGenerationError, &"did not generate the requested link after {max_tries} tries")
    var data = data0.get()
    var vertex_adjacencies = newSeq[seq[int]]()
    for (vertex, adjacencies) in data:
        if rand(0..1) == 1:
            vertex_adjacencies.add(concat(adjacencies[1 ..< adjacencies.len], adjacencies[0 ..< 1]))
        else:
            vertex_adjacencies.add(adjacencies)
    # echo vertex_adjacencies
    var edges_u = newSeq[(int, int)](num_crossings*2)
    var edges_v = newSeq[(int, int)](num_crossings*2)
    for (v, edges) in enumerate(vertex_adjacencies):
        for (i, edge) in enumerate(edges):
            if edge > 0:
                edges_u[edge-1] = (v, i)
            else:
                edges_v[-edge-1] = (v, i)
    # echo edges_u
    # echo edges_v
    var pd_code = newSeq[Crossing](num_crossings)
    for (k, strand_pair) in enumerate(zip(edges_u, edges_v)):
        let (strand1, strand2) = strand_pair
        let (a, i) = strand1
        let (b, j) = strand2
        pd_code[a][i] = k
        pd_code[b][j] = k
    return link_from_PD_code_with_extra_info[T](pd_code)
