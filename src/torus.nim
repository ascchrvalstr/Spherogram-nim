import std/strformat
import std/sequtils
import std/sugar

import links

proc torus_link_with_extra_info*[T](p0: int, q0: int): Link[T] =
    discard """
    Return a `(p,q)`-torus link, as an instance of the ``Link`` class.

    If exactly one of `p` and `q` is negative, it returns the mirror of
    `T(|p|, |q|)`.
    """
    if p0 == 0 or q0 == 0:
        raise newException(ValueError, &"torus link (p, q) requires p and q non-zero, but p = {p0} and q = {q0}")
    let p = abs(p0)
    let q = abs(q0)
    if p == 1 or q == 1:
        return link_from_PD_code(@[], name = &"T({p0},{q0})")
    let to_mirror = p0*q0 < 0
    var cur_edge_label = 0
    if p == 2:
        var pd_code = newSeq[Crossing](q)
        # set up conditions true for all two strand situations if q > 1
        pd_code[0][0] = cur_edge_label
        pd_code[q-1][1] = cur_edge_label
        cur_edge_label += 1
        pd_code[0][3] = cur_edge_label
        pd_code[q-1][2] = cur_edge_label
        cur_edge_label += 1
        pd_code[0][1] = cur_edge_label
        pd_code[1][0] = cur_edge_label
        cur_edge_label += 1
        pd_code[0][2] = cur_edge_label
        pd_code[1][3] = cur_edge_label
        cur_edge_label += 1
        # set up in between crossings
        for i in 1 ..< q-1:
            pd_code[i][1] = cur_edge_label
            pd_code[i+1][0] = cur_edge_label
            cur_edge_label += 1
            pd_code[i][2] = cur_edge_label
            pd_code[i+1][3] = cur_edge_label
            cur_edge_label += 1
        var link = link_from_PD_code(pd_code)
        if to_mirror:
            link = link.mirror()
        link.name = &"T({p0},{q0})"
        return link

    var pd_code = newSeqWith(q, newSeq[Crossing](p-1))

    # set up connecting ends
    pd_code[0][0][3] = cur_edge_label
    pd_code[q-1][0][2] = cur_edge_label
    cur_edge_label += 1
    pd_code[0][p-2][0] = cur_edge_label
    pd_code[q-1][p-2][1] = cur_edge_label
    cur_edge_label += 1

    # middle strands of connecting ends
    for i in 0 ..< p-2:
        pd_code[0][i][0] = cur_edge_label
        pd_code[q-1][i+1][2] = cur_edge_label
        cur_edge_label += 1

    if q > 1:
        # set up side connections
        for i in 0 ..< q-1:
            pd_code[i][0][2] = cur_edge_label
            pd_code[i+1][0][3] = cur_edge_label
            cur_edge_label += 1
            pd_code[i][p-2][1] = cur_edge_label
            pd_code[i+1][p-2][0] = cur_edge_label
            cur_edge_label += 1

        # set up connections between crossings
        for i in 0 ..< p-2:
            for j in 0 ..< q:
                pd_code[j][i][1] = cur_edge_label
                pd_code[j][i+1][3] = cur_edge_label
                cur_edge_label += 1

        for i in 1 ..< p-1:
            for j in 0 ..< q-1:
                pd_code[j][i][2] = cur_edge_label
                pd_code[j+1][i-1][0] = cur_edge_label
                cur_edge_label += 1

    let flattened_pd_code = collect:
        for cs in pd_code:
            for c in cs:
                c
    var link = link_from_PD_code(flattened_pd_code)
    if to_mirror:
        link = link.mirror()
    link.name = &"T({p0},{q0})"
    return link

proc torus_link*(p0: int, q0: int): Link0 =
    return torus_link_with_extra_info[int](p0, q0)
