discard """
Computing Seifert matrices and putting Links in braid closure form.

Code written by Malik Obeidin.
"""
import std/strformat
import std/sugar
import std/enumerate

import links

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
