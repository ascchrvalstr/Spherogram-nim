discard """
Computing Seifert matrices and putting Links in braid closure form.

Code written by Malik Obeidin.
"""
import std/strformat
import std/sugar

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
