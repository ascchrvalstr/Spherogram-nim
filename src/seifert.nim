discard """
Computing Seifert matrices and putting Links in braid closure form.

Code written by Malik Obeidin.
"""
import std/strformat
import std/sugar
import std/enumerate
import std/sets
import std/sequtils
import std/tables
import std/algorithm
import std/math

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

proc straighten_arrows*(arrows: var seq[array[0..3, int]]): void =
    var totally_straightened = false
    while not totally_straightened:
        totally_straightened = true
        for arrow in arrows:
            var (tail, head) = (arrow[0], arrow[1])
            if tail < head:
                # need to move tail down
                var diff = head - tail
                for (other_arrow_index, x) in enumerate(arrows):
                    if x[2] == arrow[2] and x[0] >= tail:
                        arrows[other_arrow_index][0] += diff
                    elif x[2] == arrow[2]-1 and x[1] >= tail:
                        arrows[other_arrow_index][1] += diff
                totally_straightened = false
            elif head < tail:
                # need to move head down
                var diff = tail - head
                for (other_arrow_index, x) in enumerate(arrows):
                    if x[2] == arrow[2] and x[1] >= head:
                        arrows[other_arrow_index][1] += diff
                    elif x[2] == arrow[2]+1 and x[0] >= head:
                        arrows[other_arrow_index][0] += diff
                totally_straightened = false

proc braid_arrows*[T](link0: Link[T]): seq[array[0..2, int]] =
    discard """
    Helper function to determine positions of all the crossings in a braid
    description of the link.
    """
    var link = link0.copy()
    isotope_to_braid(link)
    var circles = seifert_circles(link)
    var tree = seifert_tree(link)
    # echo circles
    # echo tree
    var tails = collect:
        for e in tree:
            e[0]
    var heads = collect(initHashSet()):
        for e in tree:
            {e[1]}
    # echo tails
    # echo heads
    var start = -1
    for (t_index, t) in enumerate(tails):
        if t notin heads:
            start = t_index
            break
    # echo start

    var ordered_strands = @[circles[start]]
    for i in 0 ..< circles.len-1:
        var new_tail = tree[start][1]
        start = tails.find(new_tail)
        ordered_strands.add(circles[start])
    # echo ordered_strands

    for i in 0 ..< ordered_strands.len-1:
        for (n, cep) in enumerate(ordered_strands[i]):
            var found_next = false
            for (m, next_cep) in enumerate(ordered_strands[i+1]):
                if cep[0] == next_cep[0]:
                    let arr = ordered_strands[i+1]
                    ordered_strands[i+1] = concat(arr[m ..< arr.len], arr[0 ..< m])
                    found_next = true
                    break
            if found_next:
                break
    # echo ordered_strands

    var positions_in_next_strand = newSeq[Table[int, (int, int)]]()
    for i in 0 ..< ordered_strands.len-1:
        var positions = initTable[int, (int, int)]()
        for (n, cep) in enumerate(ordered_strands[i]):
            for (m, next_cep) in enumerate(ordered_strands[i+1]):
                if cep[0] == next_cep[0]:
                    positions[n] = (m, (cep[1] + (if cep[1] == 2: 1 else: -1) * link.signs[cep[0]] + 4) mod 2)
                    break
        positions_in_next_strand.add(positions)
    # echo positions_in_next_strand

    var arrows = collect:
        for (n, positions) in enumerate(positions_in_next_strand):
            for (i, position) in positions.pairs:
                [i, position[0], n, position[1]]
    # echo arrows
    straighten_arrows(arrows)
    sort(arrows, proc (x, y: array[0..3, int]): int = x[0] - y[0])
    # echo arrows
    return collect:
        for arrow in arrows:
            [arrow[0], arrow[2], arrow[3]]

proc braid_word*[T](link: Link[T]): seq[int] =
    discard """
    Return a list of integers which defines a braid word whose closure is the
    given link.  The natural numbers 1, 2, 3, etc are the generators and the
    negatives are the inverses.

    Implementation follows P. Vogel, "Representation of links by
    braids, a new algorithm".
    """
    var arrows = braid_arrows(link)
    return collect:
        for arrow in arrows:
            let (_, strand, over_or_under) = (arrow[0], arrow[1], arrow[2])
            [-1, 1][over_or_under] * (strand+1)

proc seifert_matrix*[T](link: Link[T]): (seq[seq[int]], seq[seq[int]]) =
    discard """
    Returns the Seifert matrix of a link by first making it isotopic to a braid
    closure.

    Uses the algorithm described in:

    J. Collins, "An algorithm for computing the Seifert matrix of a link
    from a braid representation." (2007).
    """
    var arrows = braid_arrows(link)
    var strands = collect(initHashSet()):
        for x in arrows:
            {x[1]}
    var grouped_by_strand = collect:
        for strand in strands:
            collect:
                for x in arrows:
                    if x[1] == strand:
                        x
    var hom_gens = collect:
        for group in grouped_by_strand:
            collect:
                for i in 0 ..< group.len-1:
                    [group[i][0], group[i+1][0], group[i][2], group[i+1][2]]
    var num_gens = sum(mapIt(hom_gens, it.len))
    var matrix = newSeqWith(num_gens, newSeq[int](num_gens))
    var entries = initTable[(int, int), int]()
    var cur_index = 0
    for (i, hgi) in enumerate(hom_gens):
        for j in 0 ..< hgi.len:
            entries[(i, j)] = cur_index
            cur_index += 1
    var type_matrix = newSeqWith(num_gens, newSeq[int](num_gens))
    for (n, strand) in enumerate(hom_gens):
        # diagonal entries
        for (m, gen) in enumerate(strand):
            if gen[2] == gen[3]:
                # same sign, otherwise entry is zero
                if gen[2] == 0:
                    # both right-handed
                    matrix[entries[(n, m)]][entries[(n, m)]] = -1
                    type_matrix[entries[(n, m)]][entries[(n, m)]] = 1
                else:
                    # both left-handed
                    matrix[entries[(n, m)]][entries[(n, m)]] = 1
                    type_matrix[entries[(n, m)]][entries[(n, m)]] = 2
        # two gens on same strand, one after the other
        for (m, gen) in enumerate(strand[0 ..< strand.len-1]):
            if gen[3] == 0:
                # shared crossing is right-handed
                matrix[entries[(n, m+1)]][entries[(n, m)]] = 1
                type_matrix[entries[(n, m+1)]][entries[(n, m)]] = 3
            else:
                # shared crossing is left-handed
                matrix[entries[(n, m)]][entries[(n, m+1)]] = -1
                type_matrix[entries[(n, m)]][entries[(n, m+1)]] = 4
        # two gens on adjacent strand, "staggered"
        if n != hom_gens.len-1:
            var next_strand = hom_gens[n+1]
            for (m, gen) in enumerate(strand):
                for (l, next_gen) in enumerate(next_strand):
                    if next_gen[0] < gen[0] and gen[0] < next_gen[1] and next_gen[1] < gen[1]:
                        matrix[entries[(n+1, l)]][entries[(n, m)]] = 1
                        type_matrix[entries[(n+1, l)]][entries[(n, m)]] = 5
                    elif gen[0] < next_gen[0] and next_gen[0] < gen[1] and gen[1] < next_gen[1]:
                        matrix[entries[(n+1, l)]][entries[(n, m)]] = -1
                        type_matrix[entries[(n+1, l)]][entries[(n, m)]] = 6
    return (matrix, type_matrix)
