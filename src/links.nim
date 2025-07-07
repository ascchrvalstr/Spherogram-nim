import std/strformat
import std/tables
import std/options
import std/sugar
import std/sequtils
import std/enumerate

type Crossing* = array[0..3, int]

type LinkComponent*[T] = object of RootObj
    crossing*: int
    strand_index*: int
    extra_info*: T

proc newLinkComponent*[T](crossing: int, strand_index: int): LinkComponent[T] =
    return LinkComponent[T](crossing: crossing, strand_index: strand_index)

proc newLinkComponent*[T](crossing: int, strand_index: int, extra_info: T): LinkComponent[T] =
    return LinkComponent[T](crossing: crossing, strand_index: strand_index, extra_info: extra_info)

type Link*[T] = ref object of RootObj
    crossings*: seq[Crossing]
    signs*: seq[int]
    unlinked_unknot_components*: int
    link_components*: seq[LinkComponent[T]]
    name*: string

proc newLink*[T](crossings: seq[Crossing], signs: seq[int], unlinked_unknot_components: int, link_components: seq[LinkComponent[T]], name: string): Link[T] =
    return Link[T](crossings: crossings, signs: signs, unlinked_unknot_components: unlinked_unknot_components, link_components: link_components, name: name)

proc copy*[T](link: Link[T]): Link[T] =
    return newLink(link.crossings, link.signs, link.unlinked_unknot_components, link.link_components, link.name)

type Link0* = Link[int]

proc strand_sign(signs: seq[int], c: int, s: int): int =
    assert signs[c] in [-1, 1]
    return [-1, signs[c], 1, -signs[c]][s]

proc link_from_PD_code_with_extra_info*[T](pd_code: seq[array[0..3, int]], signs0: seq[int] = @[], unlinked_unknot_components: int = 0, name: string = ""): Link[T] =
    # pd_code and signs must have lengths equal to the number of crossings
    if signs0 != [] and pd_code.len != signs0.len:
        raise newException(ValueError, &"pd_code.len = {pd_code.len} != {signs0.len} = signs.len, where pd_code = {pd_code} and signs = {signs0}")
    # we cannot have negative number of unlinked unknot components
    if unlinked_unknot_components < 0:
        raise newException(ValueError, &"unlinked_unknot_components = {unlinked_unknot_components} < 0")
    let num_crossings = pd_code.len
    var labels = initTable[int, Option[int]]()
    var crossings = newSeq[Crossing](num_crossings)
    # glue the two strands of the same edge together, throwing an error if PD code is inconsistent
    for c in 0 ..< num_crossings:
        for s in 0 ..< 4:
            let label = pd_code[c][s]
            if label notin labels:
                labels[label] = some(c*4 + s)
            elif labels[label].isNone:
                raise newException(ValueError, &"PD code is inconsistent because label {label} appeared the third time at crossing {c}, strand {s}")
            else:
                let old_label = labels[label].get()
                crossings[c][s] = old_label
                crossings[old_label div 4][old_label mod 4] = c*4 + s
                labels[label] = none(int)
    for _, label in labels:
        if label.isSome:
            raise newException(ValueError, &"PD code is inconsistent because label {label.get()} only appeared once")
    # if signs not given, try to deduce signs
    var signs: seq[int]
    if signs0 != []:
        signs = signs0
    else:
        signs = newSeq[int](num_crossings)
        var signs_consistent = true
        var visited = newSeq[array[0..3, bool]](num_crossings)
        for c in 0 ..< num_crossings:
            if not visited[c][2]:
                var cur_c = c
                var cur_s = 2
                while not visited[cur_c][cur_s]:
                    visited[cur_c][cur_s] = true
                    visited[cur_c][(cur_s+2) mod 4] = true
                    let opposite = crossings[cur_c][cur_s]
                    # opposite should be an incoming strand
                    if opposite mod 4 == 2:
                        signs_consistent = false
                        break
                    elif opposite mod 4 in [1, 3]:
                        let proposed_sign = (if opposite mod 4 == 1: -1 else: 1)
                        if signs[opposite div 4] == -proposed_sign:
                            # the signs cannot be made consistent
                            signs_consistent = false
                            break
                        signs[opposite div 4] = proposed_sign
                    cur_c = opposite div 4
                    cur_s = (opposite+2) mod 4
                if not signs_consistent:
                    break
        # every link component must have signs assigned
        if signs_consistent:
            for sign in signs:
                if sign == 0:
                    signs_consistent = false
                    break
        if not signs_consistent:
            signs = newSeq[int](num_crossings)
    # we must either have all strands oriented, or all strands unoriented, not a mixture of the two
    var signed = -1
    for sign in signs:
        if sign notin [-1, 0, 1]:
            raise newException(ValueError, &"sign = {sign} is illegal")
        if signed == 1 and sign == 0:
            raise newException(ValueError, &"the link is known to be oriented, but we are now given sign 0")
        elif signed == 0 and sign in [-1, 1]:
            raise newException(ValueError, &"the link is known to be unoriented, but we are now given sign {sign}")
        elif sign == 0:
            signed = 0
        else:
            signed = 1
    # check that the orientations of the edges specified by the signs are consistent
    if signed == 1:
        for c in 0 ..< num_crossings:
            for s in 0 ..< 4:
                if crossings[c][s] > c*4 + s:
                    let opposite = crossings[c][s]
                    if strand_sign(signs, c, s) == strand_sign(signs, opposite div 4, opposite mod 4):
                        raise newException(ValueError, &"strand {s} of crossing {c} has the same orientation as its opposite, strand {opposite mod 4} of crossing {opposite div 4}")
    # build the link components
    var visited = newSeq[array[0..3, bool]]()
    for sign in signs:
        visited.add([sign != 0, sign == -1, false, sign == 1])
    var link_components = newSeq[LinkComponent[T]]()
    for c in 0 ..< num_crossings:
        for s in 0 ..< 4:
            if not visited[c][s]:
                link_components.add(newLinkComponent[T](c, s))
                var cur_c = c
                var cur_s = s
                while not visited[cur_c][cur_s]:
                    visited[cur_c][cur_s] = true
                    visited[cur_c][(cur_s+2) mod 4] = true
                    let opposite = crossings[cur_c][cur_s]
                    cur_c = opposite div 4
                    cur_s = (opposite+2) mod 4
    return newLink(crossings, signs, unlinked_unknot_components, link_components, name)

proc link_from_PD_code*(pd_code: seq[array[0..3, int]], signs0: seq[int] = @[], unlinked_unknot_components: int = 0, name: string = ""): Link0 =
    return link_from_PD_code_with_extra_info[int](pd_code, signs0, unlinked_unknot_components, name)

proc PD_code*[T](link: Link[T], min_strand_index: int = 0): seq[array[0..3, int]] =
    let num_crossings = link.crossings.len
    var code = newSeq[array[0..3, int]](num_crossings)
    var cur_label = min_strand_index
    var visited = newSeq[array[0..3, bool]](num_crossings)
    for comp in link.link_components:
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        while not visited[cur_c][cur_s]:
            code[cur_c][cur_s] = cur_label
            visited[cur_c][cur_s] = true
            let opposite = link.crossings[cur_c][cur_s]
            code[opposite div 4][opposite mod 4] = cur_label
            visited[opposite div 4][opposite mod 4] = true
            cur_label += 1
            cur_c = opposite div 4
            cur_s = (opposite+2) mod 4
    return code

proc gauss_code*[T](link: Link[T]): seq[seq[int]] =
    let num_crossings = link.crossings.len
    var code = newSeq[seq[int]]()
    var visited = newSeq[array[0..3, bool]](num_crossings)
    for comp in link.link_components:
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        var comp_code = newSeq[int]()
        while not visited[cur_c][cur_s]:
            let over_under = [-1, 1][cur_s mod 2]
            comp_code.add(over_under * (cur_c+1))
            visited[cur_c][cur_s] = true
            visited[cur_c][(cur_s+2) mod 4] = true
            let opposite = link.crossings[cur_c][cur_s]
            cur_c = opposite div 4
            cur_s = (opposite+2) mod 4
        code.add(comp_code)
    return code

proc extended_gauss_code*[T](link: Link[T]): seq[seq[int]] =
    let num_crossings = link.crossings.len
    var code = newSeq[seq[int]]()
    var visited = newSeq[array[0..3, bool]](num_crossings)
    var crossing_visited = newSeq[bool](num_crossings)
    for comp in link.link_components:
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        var comp_code = newSeq[int]()
        while not visited[cur_c][cur_s]:
            if not crossing_visited[cur_c]:
                let over_under = [-1, 1][cur_s mod 2]
                comp_code.add(over_under * (cur_c+1))
                crossing_visited[cur_c] = true
            else:
                comp_code.add(link.signs[cur_c] * (cur_c+1))
            visited[cur_c][cur_s] = true
            visited[cur_c][(cur_s+2) mod 4] = true
            let opposite = link.crossings[cur_c][cur_s]
            cur_c = opposite div 4
            cur_s = (opposite+2) mod 4
        code.add(comp_code)
    return code

proc mirror*[T](link: Link[T]): Link[T] =
    let num_crossings = link.crossings.len
    var new_crossings = newSeq[array[0..3, int]](num_crossings)
    var new_signs = newSeq[int](num_crossings)
    for c in 0 ..< num_crossings:
        for s in 0 ..< 4:
            let opposite = link.crossings[c][s]
            new_crossings[c][[0, 3, 2, 1][s]] = (opposite div 4)*4 + [0, 3, 2, 1][opposite mod 4]
        new_signs[c] = -link.signs[c]
    var new_link_components = newSeq[LinkComponent[T]](link.link_components.len)
    for c in 0 ..< link.link_components.len:
        let old_component = link.link_components[c]
        new_link_components[c] = newLinkComponent(old_component.crossing, [0, 3, 2, 1][old_component.strand_index], old_component.extra_info)
    return newLink(new_crossings, new_signs, link.unlinked_unknot_components, new_link_components, link.name & ".mirror")

proc `$`*[T](link: Link[T]): string =
    let name = if link.name != "": " " & link.name else: ""
    let unlinked = if link.unlinked_unknot_components != 0: &"; {link.unlinked_unknot_components} unlinked unknot comp" else: ""
    return &"<Link{name}: {link.link_components.len} comp; {link.crossings.len} cross{unlinked}>"

proc writhe*[T](link: Link[T]): int =
    if link.signs.len > 0 and link.signs[0] == 0:
        raise newException(ValueError, &"link {link} is unoriented, so its writhe is not well-defined")
    var writhe = 0
    for s in link.signs:
        writhe += s
    return writhe

type CrossingStrand*[T] = object of RootObj
    link*: Link[T]
    crossing*: int
    strand_index*: int

proc newCrossingStrand*[T](link: Link[T], crossing: int, strand_index: int): CrossingStrand[T] =
    return CrossingStrand[T](link: link, crossing: crossing, strand_index: strand_index)

proc from_pair*[T](link: Link[T], strand: (int, int)): CrossingStrand[T] = newCrossingStrand(link, strand[0], strand[1])

proc rotate*[T](strand: CrossingStrand[T], s: int = 1): CrossingStrand[T] =
    return newCrossingStrand(strand.link, strand.crossing, (strand.strand_index+s) mod 4)

proc opposite*[T](strand: CrossingStrand[T]): CrossingStrand[T] =
    let opposite = strand.link.crossings[strand.crossing][strand.strand_index]
    return newCrossingStrand(strand.link, opposite div 4, opposite mod 4)

proc next*[T](strand: CrossingStrand[T]): CrossingStrand[T] =
    let opposite = strand.link.crossings[strand.crossing][strand.strand_index]
    return newCrossingStrand(strand.link, opposite div 4, (opposite+2) mod 4)

proc next_corner*[T](strand: CrossingStrand[T]): CrossingStrand[T] =
    let opposite = strand.link.crossings[strand.crossing][strand.strand_index]
    return newCrossingStrand(strand.link, opposite div 4, (opposite+3) mod 4)

proc previous_corner*[T](strand: CrossingStrand[T]): CrossingStrand[T] =
    let opposite = strand.link.crossings[strand.crossing][(strand.strand_index+1) mod 4]
    return newCrossingStrand(strand.link, opposite div 4, opposite mod 4)

proc oriented*[T](strand: CrossingStrand[T]): CrossingStrand[T] =
    let sign = strand.link.signs[strand.crossing]
    if sign == 0:
        raise newException(ValueError, &"the link that strand {strand} points to is unoriented, so the strand cannot be oriented")
    if [-1, sign, 1, -sign][strand.strand_index] == -1:
        return strand.opposite()
    return strand

proc is_under_crossing*[T](strand: CrossingStrand[T]): bool =
    return strand.strand_index in [0, 2]

proc is_over_crossing*[T](strand: CrossingStrand[T]): bool =
    return strand.strand_index in [1, 3]

proc `$`*[T](strand: CrossingStrand[T]): string = &"<CS {strand.crossing}, {strand.strand_index}>"

proc to_pair*[T](strand: CrossingStrand[T]): (int, int) = (strand.crossing, strand.strand_index)

proc get_link_components*[T](link: Link[T]): seq[CrossingStrand[T]] =
    return collect:
        for comp in link.link_components:
            newCrossingStrand(link, comp.crossing, comp.strand_index)

proc faces*[T](link: Link[T]): seq[seq[CrossingStrand[T]]] =
    let num_crossings = link.crossings.len
    var visited = newSeq[array[0..3, bool]](num_crossings)
    var faces = newSeq[seq[CrossingStrand[T]]]()
    for c in 0 ..< num_crossings:
        for s in 0 ..< 4:
            if not visited[c][s]:
                var cur_cs = newCrossingStrand(link, c, s)
                var cur_face = newSeq[CrossingStrand[T]]()
                while not visited[cur_cs.crossing][cur_cs.strand_index]:
                    visited[cur_cs.crossing][cur_cs.strand_index] = true
                    cur_face.add(cur_cs)
                    cur_cs = cur_cs.next_corner()
                faces.add(cur_face)
    return faces

proc is_alternating*[T](link: Link[T]): bool =
    for comp in link.get_link_components():
        var cur_strand = comp
        while true:
            if cur_strand.is_over_crossing() == cur_strand.next().is_over_crossing():
                return false
            cur_strand = cur_strand.next()
            if cur_strand == comp:
                break
    return true

proc split_link_diagram*[T](link: Link[T]): seq[Link[T]] =
    let num_crossings = link.crossings.len
    # label each crossing with the index of the connected component that it belongs to
    var crossings_cc = newSeqWith(num_crossings, -1)
    # the total number of connected components
    var tot_cc = 0
    for c in 0 ..< num_crossings:
        if crossings_cc[c] == -1:
            var dfs_stack = @[c]
            crossings_cc[c] = tot_cc
            while dfs_stack.len > 0:
                let cur_crossing = dfs_stack.pop()
                for s in 0 ..< 4:
                    let adj_crossing = link.crossings[cur_crossing][s] div 4
                    if crossings_cc[adj_crossing] == -1:
                        crossings_cc[adj_crossing] = tot_cc
                        dfs_stack.add(adj_crossing)
            tot_cc += 1
    # calculate the new indices of each crossing within its connected component
    var cur_labels = newSeq[int](tot_cc)
    var indices_within_cc = newSeq[int](num_crossings)
    for c in 0 ..< num_crossings:
        indices_within_cc[c] = cur_labels[crossings_cc[c]]
        cur_labels[crossings_cc[c]] += 1
    # build the new crossings, signs and link_components
    var cc_crossings = newSeq[seq[Crossing]](tot_cc)
    var cc_signs = newSeq[seq[int]](tot_cc)
    var cc_link_components = newSeq[seq[LinkComponent[T]]](tot_cc)
    for c in 0 ..< num_crossings:
        var new_opposites: Crossing
        for s in 0 ..< 4:
            let opposite = link.crossings[c][s]
            new_opposites[s] = indices_within_cc[opposite div 4]*4 + opposite mod 4
        cc_crossings[crossings_cc[c]].add(new_opposites)
        cc_signs[crossings_cc[c]].add(link.signs[c])
    for comp in link.link_components:
        cc_link_components[crossings_cc[comp.crossing]].add(newLinkComponent[T](indices_within_cc[comp.crossing], comp.strand_index, comp.extra_info))
    # build the Link objects from the new crossings, signs and link_components
    var conn_comps = newSeq[Link[T]]()
    for cc in 0 ..< tot_cc:
        conn_comps.add(newLink(cc_crossings[cc], cc_signs[cc], 0, cc_link_components[cc], &"({link.name}).connected_component({cc})"))
    return conn_comps

proc opposite_strand*[T](link: Link[T], strand: (int, int)): (int, int) =
    let (c, s) = strand
    let opposite = link.crossings[c][s]
    return (opposite div 4, opposite mod 4)

proc next_strand*[T](link: Link[T], strand: (int, int)): (int, int) =
    let (c, s) = strand
    let opposite = link.crossings[c][s]
    return (opposite div 4, (opposite+2) mod 4)

proc previous_strand*[T](link: Link[T], strand: (int, int)): (int, int) =
    let (c, s) = strand
    let opposite = link.crossings[c][(s+2) mod 4]
    return (opposite div 4, opposite mod 4)

proc is_oriented*[T](link: Link[T]): bool =
    if link.crossings == []:
        raise newException(ValueError, &"empty link cannot be said to be oriented or unoriented")
    return link.signs[0] != 0

proc connected_sum*[T](link1: Link[T], link2: Link[T], strand1_0: (int, int) = (0, 0), strand2_0: (int, int) = (0, 0)): Link[T] =
    let new_name = &"({link1.name})#({link2.name})"
    if link2.crossings == []:
        return newLink(link1.crossings, link1.signs, link1.unlinked_unknot_components + link2.unlinked_unknot_components, link1.link_components, new_name)
    if link1.crossings == []:
        return newLink(link2.crossings, link2.signs, link1.unlinked_unknot_components + link2.unlinked_unknot_components, link2.link_components, new_name)
    let link1_oriented = link1.is_oriented()
    let link2_oriented = link2.is_oriented()
    if link1_oriented != link2_oriented:
        raise newException(ValueError, &"in connected_sum, link1.is_oriented() = {link1_oriented} != {link2_oriented} = link2.is_oriented()")
    let link_oriented = link1_oriented
    var strand1 = strand1_0
    var strand2 = strand2_0
    var visited = newSeq[array[0..3, bool]](link2.crossings.len)
    var (cur_c, cur_s) = strand2
    while not visited[cur_c][cur_s]:
        visited[cur_c][cur_s] = true
        visited[cur_c][(cur_s+2) mod 4] = true
        let opposite = link2.crossings[cur_c][cur_s]
        cur_c = opposite div 4
        cur_s = (opposite+2) mod 4
    var new_link2_components = filterIt(link2.link_components, not visited[it.crossing][it.strand_index])
    for i in 0 ..< new_link2_components.len:
        new_link2_components[i].crossing += link1.crossings.len
    if link_oriented:
        if strand_sign(link1.signs, strand1[0], strand1[1]) == 1:
            strand1 = link1.opposite_strand(strand1)
        if strand_sign(link2.signs, strand2[0], strand2[1]) == 1:
            strand2 = link2.opposite_strand(strand2)
    var strand1_op = link1.opposite_strand(strand1)
    var strand2_op = link2.opposite_strand(strand2)
    var new_crossings = link1.crossings & link2.crossings
    for c in 0 ..< link2.crossings.len:
        for s in 0 ..< 4:
            new_crossings[link1.crossings.len + c][s] += link1.crossings.len*4
    strand2[0] += link1.crossings.len
    strand2_op[0] += link1.crossings.len
    new_crossings[strand1[0]][strand1[1]] = strand2_op[0]*4 + strand2_op[1]
    new_crossings[strand1_op[0]][strand1_op[1]] = strand2[0]*4 + strand2[1]
    new_crossings[strand2[0]][strand2[1]] = strand1_op[0]*4 + strand1_op[1]
    new_crossings[strand2_op[0]][strand2_op[1]] = strand1[0]*4 + strand1[1]
    return newLink(new_crossings, link1.signs & link2.signs, link1.unlinked_unknot_components + link2.unlinked_unknot_components, link1.link_components & new_link2_components, new_name)

proc sublink*[T](link: Link[T], components: seq[bool]): Link[T] =
    if components.len != link.link_components.len:
        raise newException(ValueError, &"in sublink, components.len = {components.len} != {link.link_components.len} = link.link_components.len")
    # calculate the index of the link component for each strand, if the link component is selected
    let num_crossings = link.crossings.len
    var strand_to_comp = newSeqWith(num_crossings, [-1, -1, -1, -1])
    for comp in 0 ..< components.len:
        if components[comp]:
            let old_component = link.link_components[comp]
            var (cur_c, cur_s) = (old_component.crossing, old_component.strand_index)
            while strand_to_comp[cur_c][cur_s] == -1:
                strand_to_comp[cur_c][cur_s] = comp
                strand_to_comp[cur_c][(cur_s+2) mod 4] = comp
                (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
    # construct new crossing labels
    var cur_crossing_label = 0
    var new_crossing_labels = newSeqWith(num_crossings, -1)
    var unlinked = components
    for c in 0 ..< num_crossings:
        if strand_to_comp[c][0] != -1 and strand_to_comp[c][1] != -1:
            new_crossing_labels[c] = cur_crossing_label
            cur_crossing_label += 1
            unlinked[strand_to_comp[c][0]] = false
            unlinked[strand_to_comp[c][1]] = false
    # calculate the number of unlinked unknot components
    var unlinked_unknot_components = 0
    for c in unlinked:
        if c:
            unlinked_unknot_components += 1
    # construct the new crossings with signs
    var new_crossings = newSeqWith(cur_crossing_label, [-1, -1, -1, -1])
    var new_signs = newSeq[int](cur_crossing_label)
    for c in 0 ..< num_crossings:
        if new_crossing_labels[c] != -1:
            for s in 0 ..< 4:
                if new_crossings[new_crossing_labels[c]][s] == -1:
                    var (op_c, op_s) = link.opposite_strand((c, s))
                    while new_crossing_labels[op_c] == -1:
                        (op_c, op_s) = link.previous_strand((op_c, op_s))
                    new_crossings[new_crossing_labels[c]][s] = new_crossing_labels[op_c]*4 + op_s
                    new_crossings[new_crossing_labels[op_c]][op_s] = new_crossing_labels[c]*4 + s
            new_signs[new_crossing_labels[c]] = link.signs[c]
    # construct the new link components
    var new_link_components = newSeq[LinkComponent[T]]()
    for comp in 0 ..< components.len:
        if components[comp] and not unlinked[comp]:
            let old_component = link.link_components[comp]
            var (start_c, start_s) = (old_component.crossing, old_component.strand_index)
            while new_crossing_labels[start_c] == -1:
                (start_c, start_s) = link.previous_strand((start_c, start_s))
            new_link_components.add(newLinkComponent(new_crossing_labels[start_c], start_s, old_component.extra_info))
    return newLink(new_crossings, new_signs, unlinked_unknot_components, new_link_components, &"({link.name}).sublink({components})")

proc linking_number*[T](link: Link[T]): int =
    let num_crossings = link.crossings.len
    var strand_to_comp = newSeqWith(num_crossings, [-1, -1, -1, -1])
    for comp_index, comp in enumerate(link.link_components):
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        while strand_to_comp[cur_c][cur_s] == -1:
            strand_to_comp[cur_c][cur_s] = comp_index
            strand_to_comp[cur_c][(cur_s+2) mod 4] = comp_index
            (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
    var linking_number = 0
    for (strands, sign) in zip(strand_to_comp, link.signs):
        if strands[0] != strands[1]:
            linking_number += sign
    return linking_number div 2

proc check_crossing*[T](link: Link[T], cross: int) =
    let num_crossings = link.crossings.len
    if cross < 0 or cross >= link.crossings.len:
        raise newException(ValueError, &"crossing {cross} out of bounds [0,{num_crossings})")

proc check_crossing_strand*[T](link: Link[T], crossing_strand: (int, int)) =
    let (c, s) = crossing_strand
    link.check_crossing(c)
    if s < 0 or s >= 4:
        raise newException(ValueError, &"strand index {s} out of bounds [0, 4)")

proc is_planar*[T](link: Link[T]): bool =
    let num_crossings = link.crossings.len
    if num_crossings == 0:
        return true
    # test if the link has multiple connected components
    var crossings_visited = newSeq[bool](num_crossings)
    var dfs_stack = @[0]
    crossings_visited[0] = true
    var num_crossings_visited = 1
    while dfs_stack.len > 0:
        let cur_crossing = dfs_stack.pop()
        for s in 0 ..< 4:
            let adj_crossing = link.crossings[cur_crossing][s] div 4
            if not crossings_visited[adj_crossing]:
                crossings_visited[adj_crossing] = true
                dfs_stack.add(adj_crossing)
                num_crossings_visited += 1
    if num_crossings_visited < num_crossings:
        for conn_comp in link.split_link_diagram():
            if not is_planar(conn_comp):
                return false
        return true
    let euler = -num_crossings + len(link.faces())
    return euler == 2

proc pieces*[T](link: Link[T]): seq[seq[(int, int)]] =
    var start_pos = newSeq[(int, int)]()
    let num_crossings = link.crossings.len
    for c in 0 ..< num_crossings:
        if link.signs[c] == 0:
            start_pos.add((c, 0))
        start_pos.add((c, 2))
    # whether each strand has been visited
    var visited = newSeq[array[0..3, bool]](num_crossings)
    var over_arcs = newSeq[seq[(int, int)]]()
    for (start_c, start_s) in start_pos:
        if not visited[start_c][start_s]:
            var (cur_c, cur_s) = (start_c, start_s)
            var cur_arc = newSeq[(int, int)]()
            while true:
                cur_arc.add((cur_c, cur_s))
                visited[cur_c][cur_s] = true
                let (op_c, op_s) = link.opposite_strand((cur_c, cur_s))
                cur_arc.add((op_c, op_s))
                visited[op_c][op_s] = true
                (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
                if cur_s mod 2 == 0:
                    break
            over_arcs.add(cur_arc)
    return over_arcs

proc strand_sign*[T](link: Link[T], c: int, s: int): int =
    if link.signs[c] == 0:
        return 0
    return [-1, link.signs[c], 1, -link.signs[c]][s]
