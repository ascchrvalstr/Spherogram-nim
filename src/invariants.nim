import std/sequtils
import std/enumerate

import manu

import links

proc linking_matrix*[T](link: Link[T]): seq[seq[int]] =
    let num_crossings = link.crossings.len
    var strand_to_comp = newSeqWith(num_crossings, [-1, -1, -1, -1])
    for comp_index, comp in enumerate(link.link_components):
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        while strand_to_comp[cur_c][cur_s] == -1:
            strand_to_comp[cur_c][cur_s] = comp_index
            strand_to_comp[cur_c][(cur_s+2) mod 4] = comp_index
            (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
    let num_components = link.link_components.len
    var linking_matrix = newSeqWith(num_components, newSeqWith(num_components, 0))
    for (strands, sign) in zip(strand_to_comp, link.signs):
        if strands[0] != strands[1]:
            linking_matrix[strands[0]][strands[1]] += sign
            linking_matrix[strands[1]][strands[0]] += sign
    for i in 0 ..< num_components:
        for j in 0 ..< num_components:
            linking_matrix[i][j] = linking_matrix[i][j] div 2
    return linking_matrix

proc colorability_matrix*[T](link: Link[T]): seq[seq[int]] =
    let edges = link.pieces()
    var m = newSeqWith(link.crossings.len, newSeq[int](edges.len))
    for (ind, s) in enumerate(edges):
        for (c, i) in s:
            if i mod 2 == 1:
                m[c][ind] += 1
            else:
                m[c][ind] -= 1
    return m

proc determinant*[T](link: Link[T], method0: string = "goeritz"): int =
    # "method" is a reserved keyword in Nim, so the parameter has to be renamed to "method0"
    if method0 == "color":
        let M = link.colorability_matrix()
        let size = link.crossings.len-1
        var N = newSeqWith(size, newSeq[float64](size))
        for i in 0 ..< size:
            for j in 0 ..< size:
                N[i][j] = float64(M[i+1][j+1])
        echo matrix(N).det()
        # note that float to int conversion in Nim is similar to C++ in that the number is rounded towards zero, so we need to add 0.5 and then convert to int
        return int(abs(matrix(N).det()) + 0.5)
    if method0 == "goeritz":
        raise newException(AssertionDefect, "Not implemented")
    raise newException(AssertionDefect, "Not implemented")
