import std/sugar
import std/sequtils
import std/random
import std/enumerate
import std/strformat
import std/tables
import std/deques
import std/algorithm

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
                # if the specific crossing strand has opposite -2, this means that this strand is loose (i.e. not tied to an opposite strand)
                if link.crossings[cross][s] == -2:
                    continue
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

proc over_or_under_arcs*[T](link: Link[T], over: bool): seq[seq[(int, int)]] =
    let num_crossings = link.crossings.len
    var first_traversal_entry_points = newSeq[(int, int)]()
    var second_traversal_entry_points = newSeq[(int, int)]()
    for c in 0 ..< num_crossings:
        if link.signs[c] == 0:
            first_traversal_entry_points.add((c, 0))
        first_traversal_entry_points.add((c, 2))
        if link.signs[c] != -1:
            second_traversal_entry_points.add((c, 1))
        if link.signs[c] != 1:
            second_traversal_entry_points.add((c, 3))
    if not over:
        swap(first_traversal_entry_points, second_traversal_entry_points)
    # whether each strand has been visited
    var visited = newSeq[array[0..3, bool]](num_crossings)
    var arcs = newSeq[seq[(int, int)]]()
    for (start_c, start_s) in first_traversal_entry_points:
        if not visited[start_c][start_s]:
            var cur_arc = newSeq[(int, int)]()
            var (cur_c, cur_s) = (start_c, start_s)
            while true:
                cur_arc.add((cur_c, cur_s))
                visited[cur_c][cur_s] = true
                let (op_c, op_s) = link.opposite_strand((cur_c, cur_s))
                visited[op_c][op_s] = true
                (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
                if cur_s mod 2 == (if over: 0 else: 1):
                    break
            arcs.add(cur_arc)
    # cyclic arcs that lie totally above (for over) or totally below (for under) the rest of the diagram
    for (start_c, start_s) in second_traversal_entry_points:
        if not visited[start_c][start_s]:
            var cur_arc = newSeq[(int, int)]()
            var (cur_c, cur_s) = (start_c, start_s)
            while true:
                cur_arc.add((cur_c, cur_s))
                visited[cur_c][cur_s] = true
                let (op_c, op_s) = link.opposite_strand((cur_c, cur_s))
                visited[op_c][op_s] = true
                (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
                if (cur_c, cur_s) == (start_c, start_s):
                    break
            cur_arc.add((start_c, start_s))
            arcs.add(cur_arc)
    return arcs

proc randomize_within_lengths*[T](items: seq[seq[T]]): seq[seq[T]] =
    var max_length = -1
    for item in items:
        max_length = max(max_length, item.len)
    var by_lens = newSeq[seq[seq[T]]](max_length+1)
    for item in items:
        by_lens[item.len].add(item)
    var ans = newSeq[seq[T]]()
    for length in countdown(max_length, 0):
        by_lens[length].shuffle()
        for item_by_len in by_lens[length]:
            ans.add(item_by_len)
    return ans

proc pickup_arc*[T](link: Link[T], over: bool, arc: seq[(int, int)]): (seq[int], seq[int]) =
    if arc.len == 0:
        raise newException(ValueError, "arc is empty")
    for i in 1 ..< arc.len:
        if arc[i][1] mod 2 != (if over: 1 else: 0):
            raise newException(ValueError, &"over={over}, but arc[{i}][1] == {arc[i][1]} is inconsistent with over, where arc={arc}")
        if link.next_strand(arc[i-1]) != arc[i]:
            raise newException(ValueError, &"next strand of {arc[i-1]} is not {arc[i]}")
    var arc_is_cycle = arc.len > 1 and arc[0] == arc[arc.len-1]
    # echo arc_is_cycle
    # map from crossing to its position in arc
    var cross_to_pos = initTable[int, int]()
    for i in 1 ..< arc.len:
        cross_to_pos[arc[i][0]] = i-1
    # echo cross_to_pos
    # the number of crossings in the middle of the arc
    let num_middle_crossings = arc.len-1
    # echo num_middle_crossings
    var remove_comps = newSeq[bool](link.link_components.len)
    # update the link components whose crossing is among the to-be-removed crossings
    for (comp_index, comp) in enumerate(link.link_components):
        if comp.crossing in cross_to_pos:
            let (start_s, start_c) = (comp.crossing, comp.strand_index)
            var (cur_s, cur_c) = (start_s, start_c)
            while cur_s in cross_to_pos:
                (cur_s, cur_c) = link.previous_strand((cur_s, cur_c))
                if (cur_s, cur_c) == (start_s, start_c):
                    # we have encountered a cycle
                    break
            if (cur_s, cur_c) == (start_s, start_c):
                remove_comps[comp_index] = true
                link.unlinked_unknot_components += 1
            else:
                (link.link_components[comp_index].crossing, link.link_components[comp_index].strand_index) = (cur_s, cur_c)
    # echo link.link_components
    # echo remove_comps
    link.link_components = collect:
        for (remove_comp, comp) in zip(remove_comps, link.link_components):
            if not remove_comp:
                comp
    # echo link.link_components
    var elim = newSeq[int]()
    for i in 1 ..< arc.len:
        elim.add(arc[i][0])
    # echo elim
    # whether each middle crossing has been visited
    var crossing_visited = newSeq[bool](num_middle_crossings)
    # glue together the opposites of the other two strands of the crossings crossed by our arc, traversing along the link if necessary
    for (elim_index, elim_c) in enumerate(elim):
        if not crossing_visited[elim_index]:
            var (crossed_c, crossed_s) = link.opposite_strand(arc[elim_index])
            crossed_s = (crossed_s+3) mod 4
            # echo crossed_c
            # echo crossed_s
            var (start_c, start_s) = (crossed_c, crossed_s)
            while start_c in cross_to_pos:
                crossing_visited[cross_to_pos[start_c]] = true
                (start_c, start_s) = link.previous_strand((start_c, start_s))
                if (start_c, start_s) == (crossed_c, crossed_s):
                    break
            # echo start_c
            # echo start_s
            if (start_c, start_s) == (crossed_c, crossed_s):
                continue
            # glue the two uninvolved strands together
            var (end_c, end_s) = (crossed_c, crossed_s)
            while end_c in cross_to_pos:
                crossing_visited[cross_to_pos[end_c]] = true
                (end_c, end_s) = link.next_strand((end_c, end_s))
            end_s = (end_s+2) mod 4
            # echo end_c
            # echo end_s
            link.crossings[start_c][start_s] = end_c*4 + end_s
            link.crossings[end_c][end_s] = start_c*4 + start_s
    # echo link.crossings
    # remove the middle crossings
    var elim_updated = elim
    var (incoming_c, incoming_s) = arc[0]
    var (outgoing_c, outgoing_s) = link.opposite_strand(arc[arc.len-1])
    # echo incoming_c
    # echo incoming_s
    # echo outgoing_c
    # echo outgoing_s
    if not arc_is_cycle:
        link.crossings[incoming_c][incoming_s] = -2
        link.crossings[outgoing_c][outgoing_s] = -2
    for elim_c in elim:
        link.crossings[elim_c][0] = -1
    for i in 0 ..< num_middle_crossings:
        for j in i+1 ..< num_middle_crossings:
            update_crossing(elim_updated[j], link.crossings.len, elim_updated[i])
        if not arc_is_cycle:
            update_crossing(incoming_c, link.crossings.len, elim_updated[i])
            update_crossing(outgoing_c, link.crossings.len, elim_updated[i])
        remove_crossing(link, elim_updated[i])
    # echo elim_updated
    # echo link.crossings
    if arc_is_cycle:
        return (elim, @[])
    # build the dual graph, ignoring the incoming and outgoing strands of the given arc
    var temp_faces = newSeq[seq[(int, int)]]()
    var start_face = -1
    var end_face = -1
    var strand_to_face = newSeqWith(link.crossings.len, [-1, -1, -1, -1])
    for c in 0 ..< link.crossings.len:
        for s in 0 ..< 4:
            if (c, s) notin [(incoming_c, incoming_s), (outgoing_c, outgoing_s)] and strand_to_face[c][s] == -1:
                var temp_face = newSeq[(int, int)]()
                let start_strand = newCrossingStrand(link, c, s)
                var cur_strand = start_strand
                while strand_to_face[cur_strand.crossing][cur_strand.strand_index] == -1:
                    strand_to_face[cur_strand.crossing][cur_strand.strand_index] = temp_faces.len
                    temp_face.add((cur_strand.crossing, cur_strand.strand_index))
                    cur_strand = cur_strand.next_corner()
                    while (cur_strand.crossing, cur_strand.strand_index) in [(incoming_c, incoming_s), (outgoing_c, outgoing_s)]:
                        if cur_strand.crossing == incoming_c and cur_strand.strand_index == incoming_s mod 4:
                            start_face = temp_faces.len
                        if cur_strand.crossing == outgoing_c and cur_strand.strand_index == outgoing_s mod 4:
                            end_face = temp_faces.len
                        cur_strand.strand_index = (cur_strand.strand_index+3) mod 4
                temp_faces.add(temp_face)
    assert start_face != -1 and end_face != -1
    # echo temp_faces
    var temp_dual_graph = newSeq[seq[(int, (int, int))]](temp_faces.len)
    for (face_index, face) in enumerate(temp_faces):
        for (interface_c, interface_s) in face:
            let (op_c, op_s) = link.opposite_strand((interface_c, interface_s))
            temp_dual_graph[face_index].add((strand_to_face[op_c][op_s], (interface_c, interface_s)))
    # echo temp_dual_graph
    # use BFS to find a shortest path from start_face to end_face
    var distances = newSeqWith(temp_faces.len, -1)
    var previous_face = newSeq[int](temp_faces.len)
    var last_interface_edge = newSeq[(int, int)](temp_faces.len)
    distances[start_face] = 0
    var bfs_queue = toDeque([start_face])
    while bfs_queue.len > 0:
        let cur_face = bfs_queue.popFirst()
        for (next_face, interface_edge) in temp_dual_graph[cur_face]:
            if distances[next_face] == -1:
                distances[next_face] = distances[cur_face] + 1
                previous_face[next_face] = cur_face
                last_interface_edge[next_face] = interface_edge
                bfs_queue.addLast(next_face)
    # echo distances
    # echo previous_face
    # echo last_interface_edge
    var shortest_path_edges = newSeq[(int, int)]()
    var cur_face = end_face
    while cur_face != start_face:
        shortest_path_edges.add(last_interface_edge[cur_face])
        cur_face = previous_face[cur_face]
    shortest_path_edges.reverse()
    # echo shortest_path_edges
    # add the new intersections of the arc with the underlying diagram
    var (cur_c, cur_s) = (incoming_c, incoming_s)
    let incoming_strand_sign = link.strand_sign(incoming_c, incoming_s)
    var new_crossings = newSeq[int]()
    for (edge_c, edge_s) in shortest_path_edges:
        # the index of the opposite of strand (cur_c, cur_s)
        var last_strand_index = -1
        var new_cross_sign = -2
        var interface_edge_sign = link.strand_sign(edge_c, edge_s)
        if over:
            if interface_edge_sign == -1:
                last_strand_index = 1
                new_cross_sign = -incoming_strand_sign
            else:
                last_strand_index = 3
                new_cross_sign = incoming_strand_sign
        else:
            if incoming_strand_sign == -1:
                last_strand_index = 2
                new_cross_sign = interface_edge_sign
            else:
                last_strand_index = 0
                new_cross_sign = -interface_edge_sign
        let (op_c, op_s) = link.opposite_strand((edge_c, edge_s))
        link.crossings[cur_c][cur_s] = link.crossings.len*4 + last_strand_index
        link.crossings[edge_c][edge_s] = link.crossings.len*4 + (last_strand_index+1) mod 4
        link.crossings[op_c][op_s] = link.crossings.len*4 + (last_strand_index+3) mod 4
        var new_crossing = [-1, -1, -1, -1]
        new_crossing[last_strand_index] = cur_c*4 + cur_s
        new_crossing[(last_strand_index+1) mod 4] = edge_c*4 + edge_s
        new_crossing[(last_strand_index+3) mod 4] = op_c*4 + op_s
        (cur_c, cur_s) = (link.crossings.len, (last_strand_index+2) mod 4)
        new_crossings.add(link.crossings.len)
        link.crossings.add(new_crossing)
        link.signs.add(new_cross_sign)
    link.crossings[cur_c][cur_s] = outgoing_c*4 + outgoing_s
    link.crossings[outgoing_c][outgoing_s] = cur_c*4 + cur_s
    # echo link.crossings
    # echo link.signs
    return (elim, new_crossings)

proc pickup_arcs*[T](link: Link[T], over: bool): int =
    var link_copy = link.copy()
    var arcs = over_or_under_arcs(link, over)
    for arc in randomize_within_lengths(arcs):
        if len(arc) == 1:
            break
        let (elim, new_cross) = pickup_arc(link, over, arc)
        if elim.len > new_cross.len:
            return elim.len - new_cross.len
        assert elim.len == new_cross.len
        (link.crossings, link.signs, link.unlinked_unknot_components, link.link_components) = (link_copy.crossings, link_copy.signs, link_copy.unlinked_unknot_components, link_copy.link_components)
    return 0

proc pickup_simplify*[T](link: Link[T], type_III: int = 0): int =
    discard """
    Optimizes the overcrossings on a diagram, then the undercrossings,
    simplifying in between until the process stabilizes.
    """
    # echo link
    let init_num_crossings = link.crossings.len
    if type_III > 0:
        discard simplify_via_level_type_iii(link, type_III)
    else:
        discard basic_simplify(link)
    # echo link
    var stabilized = init_num_crossings == 0

    while not stabilized:
        let old_cross = link.crossings.len
        discard pickup_arcs(link, true)
        # echo link
        if type_III > 0:
            discard simplify_via_level_type_iii(link, type_III)
            # echo link
        
        discard pickup_arcs(link, false)
        # echo link
        if type_III > 0:
            discard simplify_via_level_type_iii(link, type_III)
            # echo link
        
        let new_cross = link.crossings.len
        stabilized = new_cross == 0 or new_cross == old_cross
    
    # TODO: connect sum with twists

    return init_num_crossings - link.crossings.len

proc simplify*[T](link: Link[T], mode: string = "basic", type_III_limit: int = 100): int =
    if mode == "basic":
        return basic_simplify(link)
    if mode == "level":
        return simplify_via_level_type_iii(link, type_III_limit)
    if mode == "pickup":
        return pickup_simplify(link)
    if mode == "global":
        return pickup_simplify(link, type_III_limit)
    raise newException(ValueError, &"mode={mode} is not among \"basic\", \"level\", \"pickup\" and \"global\"")
