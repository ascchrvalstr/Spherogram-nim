import std/sugar
import std/sequtils
import std/random
import std/enumerate
import std/tables

import links

proc remove_crossing*[T](link: Link[T], cross: int): void =
    discard """
    Deletes the given crossing. Assumes that it has already been
    disconnected from the rest of the link, so this just updates
    link.crossings and link.link_components.
    """
    let num_crossings = link.crossings.len
    if cross != num_crossings-1:
        link.crossings[cross] = link.crossings[num_crossings-1]
        link.signs[cross] = link.signs[num_crossings-1]
        # only update opposites and link_components if the last crossing is not similarly marked for deletion
        if link.crossings[cross][0] != -1:
            for s in 0 ..< 4:
                var (op_c, op_s) = link.opposite_strand((cross, s))
                if op_c == num_crossings-1:
                    op_c = cross
                link.crossings[op_c][op_s] = cross*4 + s
            for comp in 0 ..< link.link_components.len:
                if link.link_components[comp].crossing == num_crossings-1:
                    link.link_components[comp].crossing = cross
    discard link.crossings.pop()
    discard link.signs.pop()

proc update_changed*(changed: var seq[int], num_crossings: int, removed_cross: int): void =
    if removed_cross != num_crossings-1:
        for c in 0 ..< changed.len:
            if changed[c] == num_crossings-1:
                changed[c] = removed_cross

proc reidemeister_i*[T](link: Link[T], cross: int): (seq[int], seq[int]) =
    link.check_crossing(cross)
    for s in 0 ..< 4:
        if link.crossings[cross][s] == cross*4 + (s+1) mod 4:
            # we have found a possible Reidemeister I simplification
            # is the crossing the only crossing in its component?
            let only_crossing = link.crossings[cross][(s+2) mod 4] == cross*4 + (s+3) mod 4
            # echo only_crossing
            # update or delete the link component whose head is the crossing
            for comp in 0 ..< link.link_components.len:
                let old_component = link.link_components[comp]
                if old_component.crossing == cross:
                    if only_crossing:
                        link.link_components.delete(comp)
                        link.unlinked_unknot_components += 1
                    else:
                        # we go to a previous crossing that is not the current crossing
                        var (prev_c, prev_s) = link.previous_strand((old_component.crossing, old_component.strand_index))
                        if prev_c == cross:
                            # we backtrack one more step
                            (prev_c, prev_s) = link.previous_strand((prev_c, prev_s))
                        link.link_components[comp] = newLinkComponent[T](prev_c, prev_s, old_component.extra_info)
                    break
            var changed: seq[int]
            if not only_crossing:
                # glue together the opposites of the other two strands
                let (op1_c, op1_s) = link.opposite_strand((cross, (s+2) mod 4))
                let (op2_c, op2_s) = link.opposite_strand((cross, (s+3) mod 4))
                link.crossings[op1_c][op1_s] = op2_c*4 + op2_s
                link.crossings[op2_c][op2_s] = op1_c*4 + op1_s
                changed = if op1_c == op2_c: @[op1_c] else: @[op1_c, op2_c]
            else:
                changed = @[]
            # echo link.crossings
            # echo link.link_components
            # finally, remove the crossing
            let num_crossings = link.crossings.len
            remove_crossing(link, cross)
            update_changed(changed, num_crossings, cross)
            return (@[cross], changed)
    return (@[], @[])

proc update_crossing*(cross: var int, num_crossings: int, removed_cross: int): void =
    if removed_cross != num_crossings-1:
        if cross == num_crossings-1:
            cross = removed_cross

proc reidemeister_ii*[T](link: Link[T], cross: int): (seq[int], seq[int]) =
    link.check_crossing(cross)
    # find a possible Reidemeister II move
    var s = 0
    while s < 4:
        let (op0_c, op0_s) = link.opposite_strand((cross, s))
        let (op1_c, op1_s) = link.opposite_strand((cross, (s+1) mod 4))
        if op0_c == op1_c and op0_c != cross and op0_s == (op1_s+1) mod 4 and s mod 2 == op0_s mod 2:
            break
        s += 1
    # echo s
    if s == 4:
        # no possible Reidemeister II moves were found
        return (@[], @[])
    let (op_c, op1_s) = link.opposite_strand((cross, (s+1) mod 4))
    let op0_s = (op1_s+1) mod 4
    # echo op_c, op1_s, op0_s
    let num_crossings = link.crossings.len
    # degenerate case: the other two strands of op_c form a possible Reidemeister I simplification
    if link.crossings[op_c][(op0_s+1) mod 4] == op_c*4 + (op0_s+2) mod 4:
        discard reidemeister_i(link, op_c)
        var new_cross = cross
        update_crossing(new_cross, num_crossings, op_c)
        # echo new_cross
        let (_, changed) = reidemeister_i(link, new_cross)
        return (@[op_c, cross], changed)
    # degenerate case: the other two strands of cross form a possible Reidemeister I simplification
    if link.crossings[cross][(s+2) mod 4] == cross*4 + (s+3) mod 4:
        discard reidemeister_i(link, cross)
        var new_op = op_c
        update_crossing(new_op, num_crossings, cross)
        # echo new_op
        let (_, changed) = reidemeister_i(link, new_op)
        return (@[cross, op_c], changed)
    # determine whether the lower and upper components become unlinked unknot components after this simplification
    let lower_unknotted = link.crossings[cross][(s+3) mod 4] == op_c*4 + (op0_s+1) mod 4
    let upper_unknotted = link.crossings[cross][(s+2) mod 4] == op_c*4 + (op0_s+2) mod 4
    # echo lower_unknotted, upper_unknotted
    # move or remove link components based at the two to-be-removed crossings
    var remove_comp = newSeq[bool](link.link_components.len)
    for c in 0 ..< link.link_components.len:
        let old_component = link.link_components[c]
        if old_component.crossing in [cross, op_c]:
            # distinguish whether the strand lies in the lower or upper component
            var (cur_c, cur_s) = (old_component.crossing, old_component.strand_index)
            let in_upper = (cur_s-(if cur_c == cross: s else: op0_s)) mod 2 == 0
            if (if in_upper: upper_unknotted else: lower_unknotted):
                # we do not remove the component straight away; rather, we tag it for deletion later on
                remove_comp[c] = true
                link.unlinked_unknot_components += 1
            else:
                while (cur_c in [cross, op_c]):
                    (cur_c, cur_s) = link.previous_strand((cur_c, cur_s))
                link.link_components[c] = newLinkComponent[T](cur_c, cur_s, old_component.extra_info)
    # echo remove_comp
    # echo link.unlinked_unknot_components
    # delete the tagged components
    link.link_components = collect:
        for (comp, remove) in zip(link.link_components, remove_comp):
            if not remove:
                comp
    # echo link.link_components
    # glue the opposites of the two outermost strands for both the lower and the upper components, if the component is not unknotted by the Reidemeister II move
    var changed = newSeq[int]()
    if not lower_unknotted:
        let (left_c, left_s) = link.opposite_strand((cross, (s+3) mod 4))
        let (right_c, right_s) = link.opposite_strand((op_c, (op0_s+1) mod 4))
        link.crossings[left_c][left_s] = right_c*4 + right_s
        changed.add(left_c)
        link.crossings[right_c][right_s] = left_c*4 + left_s
        if right_c notin changed:
            changed.add(right_c)
    if not upper_unknotted:
        let (left_c, left_s) = link.opposite_strand((cross, (s+2) mod 4))
        let (right_c, right_s) = link.opposite_strand((op_c, (op0_s+2) mod 4))
        link.crossings[left_c][left_s] = right_c*4 + right_s
        if left_c notin changed:
            changed.add(left_c)
        link.crossings[right_c][right_s] = left_c*4 + left_s
        if right_c notin changed:
            changed.add(right_c)
    # echo link.crossings
    # echo changed
    # before removing either crossing, mark both of them for deletion
    link.crossings[cross][0] = -1
    link.crossings[op_c][0] = -1
    # remove the two crossings
    remove_crossing(link, cross)
    update_changed(changed, num_crossings, cross)
    var new_op = op_c
    update_crossing(new_op, num_crossings, cross)
    remove_crossing(link, new_op)
    update_changed(changed, num_crossings-1, new_op)
    # echo link.crossings
    # echo changed
    return (@[cross, op_c], changed)

proc update_to_visit*(to_visit: var seq[int], cross_pos: var seq[int], elim: var seq[int], changed: seq[int]): void =
    for i in 0 ..< elim.len:
        let cur_c = elim[i]
        if cur_c != cross_pos.len-1:
            let this_pos = cross_pos[cur_c]
            if this_pos != -1:
                # if the eliminated crossing is still in the stack, replace it by an invalid value
                to_visit[this_pos] = cross_pos.len-1
                cross_pos[cur_c] = -1
            let other_pos = cross_pos[cross_pos.len-1]
            if other_pos != -1:
                to_visit[other_pos] = cur_c
                cross_pos[cur_c] = other_pos
        # since elim refers to the original crossing indices, update the indices of the subsequent eliminated crossings
        for j in i+1 ..< elim.len:
            update_crossing(elim[j], cross_pos.len, cur_c)
        # since we have assigned the highest index to the eliminated crossing, we only need to remove the last element of cross_pos
        discard cross_pos.pop()
    for c in changed:
        if cross_pos[c] == -1:
            cross_pos[c] = to_visit.len
            to_visit.add(c)

proc basic_simplify*[T](link: Link[T]): int =
    let num_crossings = link.crossings.len
    var to_visit = collect:
        for c in 0 ..< num_crossings:
            num_crossings-1-c
    var cross_pos = collect:
        for c in 0 ..< num_crossings:
            num_crossings-1-c
    var total_removed = 0
    # echo to_visit
    # echo cross_pos
    while to_visit.len > 0:
        let cur_c = to_visit.pop()
        # echo cur_c
        if cur_c >= link.crossings.len:
            continue
        cross_pos[cur_c] = -1
        # first try to perform a Reidemeister I move
        var (elim, changed) = reidemeister_i(link, cur_c)
        # echo elim
        # echo changed
        if elim.len > 0:
            total_removed += elim.len
            update_to_visit(to_visit, cross_pos, elim, changed)
        else:
            # then try to perform a Reidemeister II move
            (elim, changed) = reidemeister_ii(link, cur_c)
            # echo elim
            # echo changed
            if elim.len > 0:
                total_removed += elim.len
                update_to_visit(to_visit, cross_pos, elim, changed)
        # echo link.crossings
        # echo to_visit
        # echo cross_pos
        # echo total_removed
    return total_removed

proc possible_type_iii_moves*[T](link: Link[T]): seq[array[0..2, (int, int)]] =
    var ans = newSeq[array[0..2, (int, int)]]()
    for face in link.faces():
        if face.len != 3:
            continue
        let cs = [face[0].crossing, face[1].crossing, face[2].crossing]
        # we must have three distinct crossings
        if cs[0] == cs[1] or cs[1] == cs[2] or cs[2] == cs[0]:
            continue
        # there must be one strand on top and another at bottom, otherwise no strand can be moved from one side of the remaining crossing to the other side
        let ss = [face[0].strand_index, face[1].strand_index, face[2].strand_index]
        if (ss[0] mod 2 + ss[1] mod 2 + ss[2] mod 2) in [0, 3]:
            continue
        # normalize the face so that the strand going through the second and third crossings is the bottom strand
        var s = 0
        while not (face[(s+1) mod 3].strand_index mod 2 == 0 and face[(s+2) mod 3].strand_index mod 2 == 1):
            s += 1
        ans.add([(cs[s], ss[s]), (cs[(s+1) mod 3], ss[(s+1) mod 3]), (cs[(s+2) mod 3], ss[(s+2) mod 3])])
    return ans

proc reidemeister_iii*[T](link: Link[T], triple: array[0..2, (int, int)]): void =
    # echo link.crossings
    # echo triple
    var cs = [triple[0][0], triple[1][0], triple[2][0]]
    var old_crossings = [link.crossings[cs[0]], link.crossings[cs[1]], link.crossings[cs[2]]]
    var ss = [triple[0][1], triple[1][1], triple[2][1]]
    # echo cs
    # echo old_crossings
    # echo ss
    link.crossings[cs[0]][ss[0]] = (
        if old_crossings[1][(ss[1]+3) mod 4] div 4 == cs[2]:
            cs[0]*4 + (ss[0]+1) mod 4
        elif old_crossings[1][(ss[1]+3) mod 4] div 4 == cs[1]:
            cs[2]*4 + (ss[2]+1) mod 4
        else:
            old_crossings[1][(ss[1]+3) mod 4]
    )
    link.crossings[cs[0]][(ss[0]+1) mod 4] = (
        if old_crossings[2][(ss[2]+2) mod 4] div 4 == cs[2]:
            cs[1]*4 + ss[1]
        elif old_crossings[2][(ss[2]+2) mod 4] div 4 == cs[1]:
            cs[0]*4 + ss[0]
        else:
            old_crossings[2][(ss[2]+2) mod 4]
    )
    link.crossings[cs[0]][(ss[0]+2) mod 4] = cs[1]*4 + (ss[1]+3) mod 4
    link.crossings[cs[0]][(ss[0]+3) mod 4] = cs[2]*4 + (ss[2]+2) mod 4
    link.crossings[cs[1]][ss[1]] = (
        if old_crossings[2][(ss[2]+3) mod 4] div 4 == cs[0]:
            cs[1]*4 + (ss[1]+1) mod 4
        elif old_crossings[2][(ss[2]+3) mod 4] div 4 == cs[2]:
            cs[0]*4 + (ss[0]+1) mod 4
        else:
            old_crossings[2][(ss[2]+3) mod 4]
    )
    link.crossings[cs[1]][(ss[1]+1) mod 4] = (
        if old_crossings[0][(ss[0]+2) mod 4] div 4 == cs[0]:
            cs[2]*4 + ss[2]
        elif old_crossings[0][(ss[0]+2) mod 4] div 4 == cs[2]:
            cs[1]*4 + ss[1]
        else:
            old_crossings[0][(ss[0]+2) mod 4]
    )
    link.crossings[cs[1]][(ss[1]+2) mod 4] = cs[2]*4 + (ss[2]+3) mod 4
    link.crossings[cs[1]][(ss[1]+3) mod 4] = cs[0]*4 + (ss[0]+2) mod 4
    link.crossings[cs[2]][ss[2]] = (
        if old_crossings[0][(ss[0]+3) mod 4] div 4 == cs[1]:
            cs[2]*4 + (ss[2]+1) mod 4
        elif old_crossings[0][(ss[0]+3) mod 4] div 4 == cs[0]:
            cs[1]*4 + (ss[1]+1) mod 4
        else:
            old_crossings[0][(ss[0]+3) mod 4]
    )
    link.crossings[cs[2]][(ss[2]+1) mod 4] = (
        if old_crossings[1][(ss[1]+2) mod 4] div 4 == cs[1]:
            cs[0]*4 + ss[0]
        elif old_crossings[1][(ss[1]+2) mod 4] div 4 == cs[0]:
            cs[2]*4 + ss[2]
        else:
            old_crossings[1][(ss[1]+2) mod 4]
    )
    link.crossings[cs[2]][(ss[2]+2) mod 4] = cs[0]*4 + (ss[0]+3) mod 4
    link.crossings[cs[2]][(ss[2]+3) mod 4] = cs[1]*4 + (ss[1]+2) mod 4
    # echo link.crossings
    for i in 0 ..< 3:
        for s in 0 ..< 4:
            let opposite = link.crossings[cs[i]][s]
            link.crossings[opposite div 4][opposite mod 4] = cs[i]*4 + s
    # echo link.crossings

proc simplify_via_level_type_iii*[T](link: Link[T], max_consecutive_failures: int = 100): int =
    discard """
    Applies a series of type III moves to the link, simplifying it via type
    I and II moves whenever possible.
    """
    var total_removed = 0
    total_removed += basic_simplify(link)
    var failures = 0
    while failures < max_consecutive_failures:
        let poss_moves = possible_type_iii_moves(link)
        if poss_moves.len == 0:
            break
        reidemeister_iii(link, sample(poss_moves))
        # echo link.crossings
        let removed = basic_simplify(link)
        # echo removed
        # echo link.crossings
        if removed > 0:
            failures = 0
            total_removed += removed
        else:
            failures += 1
    return total_removed

proc simplify*[T](link: Link[T], mode: string = "basic", type_III_limit: int = 100): int =
    if mode == "basic":
        return basic_simplify(link)
    if mode == "level":
        return simplify_via_level_type_iii(link, type_III_limit)
    raise newException(AssertionDefect, "Not implemented")

type DualGraphOfFaces* = seq[seq[(int, (int, int))]]

proc dual_graph*[T](link: Link[T]): DualGraphOfFaces =
    let faces = link.faces()
    var edge_to_face = newSeq[array[0..3, int]](link.crossings.len)
    for (face_index, face) in enumerate(faces):
        for strand in face:
            edge_to_face[strand.crossing][strand.strand_index] = face_index
    var dual = newSeq[seq[(int, (int, int))]](faces.len)
    for (face_index, face) in enumerate(faces):
        for strand in face:
            let opposite_strand = strand.opposite()
            dual[face_index].add((edge_to_face[opposite_strand.crossing][opposite_strand.strand_index], (strand.crossing, strand.strand_index)))
    return dual

proc deconnect_sum*[T](link0: Link[T]): seq[Link[T]] =
    var link = link0.copy()
    # echo link.faces()
    var dual = link.dual_graph()
    # echo dual
    let num_faces = dual.len
    for cur_face in 0 ..< num_faces-1:
        # maps destination face to the list of dual edges from this face to the destination face
        var dual_edges_by_face = initTable[int, seq[(int, int)]]()
        for (dest_face, edge) in dual[cur_face]:
            if dest_face > cur_face:
                if dest_face notin dual_edges_by_face:
                    dual_edges_by_face[dest_face] = newSeq[(int, int)]()
                dual_edges_by_face[dest_face].add(edge)
        for dest_face, edges in dual_edges_by_face.pairs:
            if edges.len > 1:
                # echo cur_face
                # echo dest_face
                # echo edges
                # figure out the common link component
                let (this_c, this_s) = edges[0]
                let (op_c, op_s) = link0.opposite_strand((this_c, this_s))
                var comp_index = 0
                var found = false
                var (cur_c, cur_s) = (-1, -1)
                while comp_index < link.link_components.len:
                    let cur_comp = link.link_components[comp_index]
                    (cur_c, cur_s) = (cur_comp.crossing, cur_comp.strand_index)
                    while true:
                        if (cur_c, cur_s) in [(this_c, this_s), (op_c, op_s)]:
                            found = true
                            break
                        (cur_c, cur_s) = link.previous_strand((cur_c, cur_s))
                        if (cur_c, cur_s) == (cur_comp.crossing, cur_comp.strand_index):
                            break
                    if found:
                        break
                    comp_index += 1
                assert found and comp_index < link.link_components.len
                # echo link.link_components
                # echo comp_index
                let extra_info = link.link_components[comp_index].extra_info
                var new_link_components = link.link_components[0 ..< comp_index]
                for i in 0 ..< edges.len:
                    var (add_c, add_s) = edges[i]
                    if (cur_c, cur_s) == (op_c, op_s):
                        (add_c, add_s) = link0.opposite_strand(edges[i])
                    new_link_components.add(newLinkComponent(add_c, add_s, extra_info))
                link.link_components = concat(new_link_components, link.link_components[comp_index+1 ..< link.link_components.len])
                # echo link.link_components
                for i in 0 ..< edges.len:
                    let (prev_c, prev_s) = link0.opposite_strand(edges[i])
                    let (next_c, next_s) = edges[(i+1) mod edges.len]
                    link.crossings[prev_c][prev_s] = next_c*4 + next_s
                    link.crossings[next_c][next_s] = prev_c*4 + prev_s
                # echo link.crossings
    return link.split_link_diagram()
