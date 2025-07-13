discard """
Switching crossings to make each twist region consistent.
"""
import ../src/links

proc make_twist_regions_consistent*[T](link: Link[T]): int =
    discard """
    Changes crossings so that no bigon permits a Type II move cancelling
    the pair of crossings at the end of the bigon.  The system is that the
    end crossing with the lowest label "wins" in determining the type of the twist
    region.  The code assumes that no link component is effectively a meridian
    loop for another component; said differently, no two bigons share a
    common edge.

    Note, this implementation fails if given something like a (2, n) torus link.
    """
    let num_crossings = link.crossings.len
    var total_switched = 0
    for cur_c in 0 ..< num_crossings:
        var s = 0
        while s < 4:
            let (op0_c, op0_s) = link.opposite_strand((cur_c, s))
            let (op1_c, op1_s) = link.opposite_strand((cur_c, (s+1) mod 4))
            if op0_c == op1_c and op0_c != cur_c and op0_s == (op1_s+1) mod 4 and s mod 2 == op0_s mod 2 and op0_c > cur_c:
                link.rotate_crossing(op0_c, [3, 1, 1][link.signs[op0_c]+1])
                total_switched += 1
                s += 1
            s += 1
    return total_switched
