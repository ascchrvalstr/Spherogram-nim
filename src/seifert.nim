discard """
Computing Seifert matrices and putting Links in braid closure form.

Code written by Malik Obeidin.
"""
import std/strformat
import std/sugar
import std/enumerate
import std/sets
import std/sequtils

import links
import simplify

proc seifert_circles*[T](link: Link[T]): seq[seq[(int, int)]] =
    discard """
    Returns the circles in the diagram created by Seifert's algorithm
    """
    if link.crossings.len > 0 and link.signs[0] == 0:
        raise newException(ValueError, &"link {link} is unoriented, but Seifert's algorithm is only defined for oriented links")
    let num_crossings = link.crossings.len
    var visited = collect:
        for sign in link.signs:
            [true, sign == -1, false, sign == 1]
    var circles = newSeq[seq[(int, int)]]()
    for c in 0 ..< num_crossings:
        for s in 0 ..< 4:
            if not visited[c][s]:
                var circle = newSeq[(int, int)]()
                var (cur_c, cur_s) = (c, s)
                while not visited[cur_c][cur_s]:
                    visited[cur_c][cur_s] = true
                    circle.add((cur_c, cur_s))
                    let (op_c, op_s) = link.opposite_strand((cur_c, cur_s))
                    cur_c = op_c
                    if op_s == 0:
                        cur_s = (op_s + link.signs[op_c] + 4) mod 4
                    else:
                        cur_s = (op_s - link.signs[op_c] + 4) mod 4
                circles.add(circle)
    return circles

proc admissible_moves*[T](link: Link[T]): (seq[((int, int), (int, int))], seq[(int, int)]) =
    var circles = seifert_circles(link)
    let num_crossings = link.crossings.len
    var cs_to_seifert_circle = newSeq[array[0..3, int]](num_crossings)
    for (circle_index, circle) in enumerate(circles):
        for (strand_c, strand_s) in circle:
            cs_to_seifert_circle[strand_c][strand_s] = circle_index
            let (op_c, op_s) = link.opposite_strand((strand_c, strand_s))
            cs_to_seifert_circle[op_c][op_s] = circle_index
    # echo cs_to_seifert_circle
    var pairs = newSeq[((int, int), (int, int))]()
    var seifert_circle_pairs = newSeq[(int, int)]()
    for face in link.faces():
        for i in 0 ..< face.len-1:
            for j in i+1 ..< face.len:
                let (cs1, cs2) = (face[i], face[j])
                let (circle1, circle2) = (cs_to_seifert_circle[cs1.crossing][cs1.strand_index],
                                          cs_to_seifert_circle[cs2.crossing][cs2.strand_index])
                if circle1 != circle2:
                    pairs.add((cs1.to_pair(), cs2.to_pair()))
                    seifert_circle_pairs.add((circle1, circle2))
    return (pairs, seifert_circle_pairs)

proc seifert_tree*[T](link: Link[T]): seq[array[0..1, HashSet[int]]] =
    discard """
    The oriented tree corresponding to the complementary regions of
    the Seifert circles.
    """
    var circles = seifert_circles(link)
    # echo circles
    let num_crossings = link.crossings.len
    var strand_to_circle = newSeqWith(num_crossings, [-1, -1, -1, -1])
    for (circle_index, circle) in enumerate(circles):
        for (strand_c, strand_s) in circle:
            strand_to_circle[strand_c][strand_s] = circle_index
    # echo strand_to_circle
    var edges = collect:
        for n in 0 ..< circles.len:
            [toHashSet([n]), toHashSet([n])]
    # echo edges
    for c in 0 ..< num_crossings:
        var under_circle = strand_to_circle[c][if link.signs[c] == 1: 1 else: 3]
        var over_circle = strand_to_circle[c][2]
        # echo under_circle, " ", over_circle
        if link.signs[c] == -1:
            swap(under_circle, over_circle)
        edges[over_circle][1] = union(edges[over_circle][1], edges[under_circle][0])
        edges[under_circle][0] = edges[over_circle][1]
        # echo edges
    # connect all vertices which intersect
    let num_circles = circles.len
    for e1_index in 0 ..< num_circles-1:
        let e1 = edges[e1_index]
        for e2_index in e1_index+1 ..< num_circles:
            let e2 = edges[e2_index]
            for i in 0 ..< 2:
                for j in 0 ..< 2:
                    if intersection(e1[i], e2[j]).len > 1:
                        edges[e1_index][i] = union(edges[e1_index][i], edges[e2_index][j])
                        edges[e2_index][j] = edges[e1_index][i]
    return edges

proc remove_admissible_move*[T](link: Link[T]): bool =
    discard """
    Performs a Reidemester II move to remove one branching point of the Seifert
    tree. The goal is to turn the Seifert tree into a chain.
    """
    var (moves, circle_pairs) = admissible_moves(link)
    var tree = seifert_tree(link)
    # echo moves
    # echo circle_pairs
    # echo tree
    for e1_index in 0 ..< tree.len-1:
        let e1 = tree[e1_index]
        for e2_index in e1_index+1 ..< tree.len:
            let e2 = tree[e2_index]
            if e1[0] == e2[0] or e1[1] == e2[1]:
                # edges start or end at the same point
                # echo e1_index, " ", e2_index
                for (n, pair) in enumerate(circle_pairs):
                    if toHashSet([pair[0], pair[1]]) == toHashSet([e1_index, e2_index]):
                        # echo n, " ", pair
                        var (cs1, cs2) = moves[n]
                        # echo cs1, " ", cs2
                        reverse_type_ii(link, cs2, cs1)
                        return true
    return false

proc isotope_to_braid*[T](link: Link[T]): void =
    discard """
    Performs Reidemester II moves until the Seifert tree becomes a chain, i.e.
    the Seifert circles are all nested and compatibly oriented, following
    P. Vogel, "Representation of links by braids, a new algorithm"
    """
    while remove_admissible_move(link):
        discard

proc is_chain*(tree: seq[array[0..1, HashSet[int]]]): bool =
    var tails = collect:
        for e in tree:
            e[0]
    var heads = collect:
        for e in tree:
            e[1]
    return toHashSet(tails).len == tails.len and toHashSet(heads).len == heads.len
