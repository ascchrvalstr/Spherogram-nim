import std/unittest

import ../src/links
import ../src/simplify

test "update_changed":
    # removed_crossing == num_crossings-1
    var changed = @[0, 1]
    update_changed(changed, 3, 2)
    check changed == [0, 1]
    # no last crossing in changed
    changed = @[0, 2]
    update_changed(changed, 4, 1)
    check changed == [0, 2]
    # has last crossing in changed
    changed = @[3, 0]
    update_changed(changed, 4, 1)
    check changed == [1, 0]

test "reidemeister_i":
    # check that the one vertex unknot is correctly simplified
    let one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    check reidemeister_i(one_vertex_unknot, 0) == (@[0], @[])
    # check only_crossing == true
    check one_vertex_unknot.crossings == []
    check one_vertex_unknot.signs == []
    check one_vertex_unknot.unlinked_unknot_components == 1
    check one_vertex_unknot.link_components == []

    # check that remove_crossing is correctly implemented, especially regarding relabeling of gluings
    let two_unknots = link_from_PD_code(@[[0, 1, 1, 0], [2, 2, 3, 3]])
    check reidemeister_i(two_unknots, 0) == (@[0], @[])
    # check link.crossings == @[[3, 2, 1, 0], [5, 4, 7, 6]]
    # check link.link_components == @[newLinkComponent[int](1, 1)]
    check two_unknots.crossings == [[1, 0, 3, 2]]
    check two_unknots.signs == [1]
    check two_unknots.unlinked_unknot_components == 1
    check two_unknots.link_components == [newLinkComponent[int](0, 1)]

    # check that the opposites of the other two strands are correctly glued
    let two_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 2], [2, 3, 3, 0]])
    check reidemeister_i(two_vertex_unknot, 0) == (@[0], @[0])
    # check only_crossing == false
    # check link.crossings == @[[7, 2, 1, 4], [7, 6, 5, 4]]
    # check link.link_components == @[newLinkComponent[int](1, 3)]
    check two_vertex_unknot.crossings == [[3, 2, 1, 0]]
    check two_vertex_unknot.signs == [-1]
    check two_vertex_unknot.unlinked_unknot_components == 0
    check two_vertex_unknot.link_components == [newLinkComponent[int](0, 3)]

    # check that we are backtracking once more if the link component still starts with the crossing to be removed after one round of backtracking
    let two_vertex_unknot_2 = link_from_PD_code(@[[0, 1, 1, 2], [2, 3, 3, 0]])
    two_vertex_unknot_2.link_components = @[newLinkComponent[int](0, 3)]
    check reidemeister_i(two_vertex_unknot_2, 0) == (@[0], @[0])
    # check link.link_components == @[newLinkComponent[int](1, 3)]
    check two_vertex_unknot_2.link_components == [newLinkComponent[int](0, 3)]

    # check that (@[], @[]) is returned if no possible Reidemeister I moves are found
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    for c in 0 ..< 3:
        check reidemeister_i(trefoil, c) == (@[], @[])

test "update_crossing":
    # removed_crossing == num_crossings-1
    var cross = 1
    update_crossing(cross, 3, 2)
    check cross == 1
    # not last crossing
    cross = 0
    update_crossing(cross, 3, 1)
    check cross == 0
    # is last crossing
    cross = 2
    update_crossing(cross, 3, 1)
    check cross == 1

test "reidemeister_ii":
    # check that (@[], @[]) is returned if no possible Reidemeister II moves are found
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    for c in 0 ..< 3:
        check reidemeister_ii(trefoil, c) == (@[], @[])
        # check s == 4

    # check the right degenerate case
    var right_degenerate = link_from_PD_code(@[[5, 3, 0, 2], [1, 1, 2, 0], [4, 3, 5, 4]])
    check reidemeister_ii(right_degenerate, 0) == (@[2, 0], @[0])
    # check s == 0
    # check op_c == 2 and op1_s == 1 and op2_s == 2
    # check new_cross == 0
    check right_degenerate.crossings == [[1, 0, 3, 2]]
    check right_degenerate.signs == [1]
    check right_degenerate.unlinked_unknot_components == 0
    check right_degenerate.link_components == [newLinkComponent[int](0, 2)]
    # check that in the right degenerate case, new_cross is updated if cross is moved by the removal of op_c
    right_degenerate = link_from_PD_code(@[[4, 3, 5, 4], [1, 1, 2, 0], [5, 3, 0, 2]])
    check reidemeister_ii(right_degenerate, 2) == (@[0, 2], @[0])
    # check op_c == 0
    # check new_cross == 0

    # check the left degenerate case
    let left_degenerate = link_from_PD_code(@[[5, 1, 0, 0], [4, 1, 5, 2], [2, 3, 3, 4]])
    check reidemeister_ii(left_degenerate, 0) == (@[0, 1], @[0])
    # check new_op == 1

    # check the simplest unknotted link with two components
    var unknotted_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 3, 1, 2]], @[-1, 1])
    check reidemeister_ii(unknotted_hopf, 0) == (@[0, 1], @[])
    # check lower_unknotted and upper_unknotted
    # check remove_comp == @[true, true]
    # check link.link_components == @[]
    # check link.crossings == @[[6, 5, 4, 7], [2, 1, 0, 3]]
    # check changed == @[]
    # check link.crossings == @[]
    # check changed == @[]
    unknotted_hopf = link_from_PD_code(@[[0, 4, 1, 5], [3, 4, 0, 5], [1, 3, 2, 2]], @[-1, 1, 1])
    check reidemeister_ii(unknotted_hopf, 0) == (@[0, 1], @[0])
    # check lower_unknotted and not upper_unknotted
    # check remove_comp == @[false, true]
    # check link.link_components == @[newLinkComponent[int](2, 1)]
    # check link.crossings == @[[6, 5, 8, 7], [9, 1, 0, 3], [9, 8, 11, 10]]
    # check changed == @[2]
    # check link.crossings == @[[1, 0, 3, 2]]
    # check changed == @[0]
    unknotted_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 5, 1, 2], [3, 4, 4, 5]])
    check reidemeister_ii(unknotted_hopf, 0) == (@[0, 1], @[0])
    # check not lower_unknotted and upper_unknotted
    # check remove_comp == @[false, true]
    # check link.link_components == @[newLinkComponent[int](2, 3)]
    # check link.crossings == @[[6, 8, 4, 7], [2, 11, 0, 3], [11, 10, 9, 8]]
    # check changed == @[2]
    # check link.crossings == @[[3, 2, 1, 0]]
    # check changed == @[0]

    # check that in case the second removed crossing is the last one, so that it takes the place of the first removed crossing, its opposite strands are not updated because the second removed crossing has been marked for deletion
    unknotted_hopf = link_from_PD_code(@[[0, 1, 1, 2], [2, 5, 3, 4], [5, 6, 6, 7], [3, 7, 0, 4]])
    check reidemeister_ii(unknotted_hopf, 1) == (@[1, 3], @[1, 0])
    check unknotted_hopf.crossings == [[3, 2, 1, 0], [7, 6, 5, 4]]

test "basic_simplify":
    let trefoil_changed_crossing = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [4, 2, 5, 1]])
    check basic_simplify(trefoil_changed_crossing) == 3
    check trefoil_changed_crossing.crossings == []
    check trefoil_changed_crossing.signs == []
    check trefoil_changed_crossing.unlinked_unknot_components == 1
    check trefoil_changed_crossing.link_components == []

test "possible_type_iii_moves":
    let degenerate_3edge_face = link_from_PD_code(@[[0, 3, 1, 0], [2, 2, 3, 1]])
    check possible_type_iii_moves(degenerate_3edge_face) == []
    let impossible_triangle = link_from_PD_code(@[[0, 1, 1, 2], [2, 3, 3, 4], [4, 5, 5, 0]])
    check possible_type_iii_moves(impossible_triangle) == []
    let unordered_type_iii = link_from_PD_code(@[[5, 5, 0, 4], [0, 1, 1, 2], [3, 3, 4, 2]])
    check possible_type_iii_moves(unordered_type_iii) == [[(2, 2), (0, 2), (1, 3)]]
    let ordered_type_iii = link_from_PD_code(@[[3, 3, 4, 2], [5, 5, 0, 4], [0, 1, 1, 2]])
    check possible_type_iii_moves(ordered_type_iii) == [[(0, 2), (1, 2), (2, 3)]]

test "reidemeister_iii":
    let same_crossing_glued = link_from_PD_code(@[[3, 3, 4, 2], [5, 5, 0, 4], [0, 1, 1, 2]])
    reidemeister_iii(same_crossing_glued, [(0, 2), (1, 2), (2, 3)])
    check same_crossing_glued.crossings == @[[5, 9, 8, 6], [10, 0, 3, 11], [2, 1, 4, 7]]
    let adjacent_crossings_glued = link_from_PD_code(@[[1, 5, 2, 4], [5, 3, 0, 2], [0, 3, 1, 4]])
    reidemeister_iii(adjacent_crossings_glued, [(0, 2), (1, 2), (2, 3)])
    check adjacent_crossings_glued.crossings == @[[5, 9, 3, 2], [10, 0, 7, 6], [11, 1, 4, 8]]
    let all_strands_different = link_from_PD_code(@[[0, 3, 1, 4], [2, 1, 3, 2], [7, 5, 8, 4], [5, 7, 6, 6], [11, 9, 0, 8], [9, 11, 10, 10]])
    reidemeister_iii(all_strands_different, [(2, 2), (4, 2), (0, 3)])
    check all_strands_different.crossings == @[[21, 9, 16, 12], [7, 18, 11, 4], [17, 1, 20, 6], [3, 19, 15, 14], [2, 8, 5, 13], [10, 0, 23, 22]]

    # new special case detected by running RL-for-KnotTheory
    # one strand is glued to its third next neighbor, ccw or cw (which yield the same strand)
    let glued_to_third_neighbor = link_from_PD_code(@[[1, 5, 0, 2], [2, 0, 3, 1], [3, 4, 4, 5]])
    reidemeister_iii(glued_to_third_neighbor, [(0, 0), (1, 2), (2, 3)])
    check glued_to_third_neighbor.crossings == @[[7, 6, 5, 9], [10, 2, 1, 0], [11, 3, 4, 8]]

test "simplify_via_level_type_iii":
    # 13n3370, with one crossing change via band attachment, with another band attached (and thus the diagram of an unknot)
    var cc2 = link_from_PD_code(@[[0, 21, 1, 22], [20, 1, 21, 2], [2, 19, 36, 20], [18, 32, 19, 13], [17, 30, 18, 27], [25, 16, 26, 17], [4, 16, 5, 15], [14, 4, 15, 3], [7, 13, 8, 12], [11, 28, 12, 29], [10, 24, 11, 23], [22, 10, 23, 9], [29, 8, 0, 9], [6, 27, 7, 28], [24, 5, 25, 6], [26, 14, 33, 31], [30, 31, 34, 32], [33, 3, 37, 35], [35, 37, 36, 34]])
    check simplify_via_level_type_iii(cc2) == 19
    check cc2.crossings == []
    check cc2.unlinked_unknot_components == 1
    # performance:
    # simplify_via_level_type_iii(cc2): ~1500 per second vs ~150 for Spherogram

test "dual_graph":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.dual_graph() == [@[(3, (0, 0)), (1, (2, 2))], @[(0, (0, 1)), (4, (2, 1)), (2, (1, 1))], @[(1, (0, 2)), (3, (1, 0))], @[(2, (0, 3)), (4, (1, 3)), (0, (2, 3))], @[(1, (1, 2)), (3, (2, 0))]]
    # check the one vertex unknot
    let one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    check one_vertex_unknot.dual_graph() == @[@[(2, (0, 0)), (1, (0, 2))], @[(0, (0, 1))], @[(0, (0, 3))]]
    # check the case where the dual graph has multiple edges between the same pair of faces
    let two_vertex_unknot = link_from_PD_code(@[[0, 3, 1, 0], [2, 1, 3, 2]])
    check two_vertex_unknot.dual_graph() == [@[(2, (0, 0)), (1, (0, 2)), (3, (1, 0)), (1, (1, 2))], @[(0, (0, 1)), (0, (1, 1))], @[(0, (0, 3))], @[(0, (1, 3))]]

test "deconnect_sum":
    let two_vertex_unknot = link_from_PD_code(@[[0, 3, 1, 0], [2, 1, 3, 2]])
    var result = two_vertex_unknot.deconnect_sum()
    # check link.link_components == @[(crossing: 0, strand_index: 2, extra_info: 0), (crossing: 1, strand_index: 2, extra_info: 0)]
    # check link.crossings == @[[3, 2, 1, 0], [7, 6, 5, 4]]
    check result.len == 2
    check result[0].crossings == [[3, 2, 1, 0]]
    check result[0].signs == [-1]
    check result[0].link_components == [newLinkComponent[int](0, 2)]
    check result[1].crossings == [[3, 2, 1, 0]]
    check result[1].signs == [-1]
    check result[1].link_components == [newLinkComponent[int](0, 2)]
    let three_vertex_unknot = link_from_PD_code(@[[0, 2, 1, 1], [4, 0, 5, 5], [2, 4, 3, 3]])
    result = three_vertex_unknot.deconnect_sum()
    # check link.link_components == @[(crossing: 1, strand_index: 1, extra_info: 0), (crossing: 2, strand_index: 1, extra_info: 0), (crossing: 0, strand_index: 1, extra_info: 0)]
    # check link.crossings == @[[1, 0, 3, 2], [5, 4, 7, 6], [9, 8, 11, 10]]
    check result.len == 3
    check result[0].crossings == [[1, 0, 3, 2]]
    check result[0].signs == [1]
    check result[0].link_components == [newLinkComponent[int](0, 1)]
    check result[1].crossings == [[1, 0, 3, 2]]
    check result[1].signs == [1]
    check result[1].link_components == [newLinkComponent[int](0, 1)]
    check result[2].crossings == [[1, 0, 3, 2]]
    check result[2].signs == [1]
    check result[2].link_components == [newLinkComponent[int](0, 1)]
    # trefoil \# one_vertex_unknot
    let trefoil_unknot = link_from_PD_code(@[[7, 2, 0, 3], [3, 0, 4, 1], [5, 5, 6, 4], [1, 6, 2, 7]])
    result = trefoil_unknot.deconnect_sum()
    # check link.link_components == @[(crossing: 2, strand_index: 2, extra_info: 0), (crossing: 1, strand_index: 2, extra_info: 0)]
    # check link.crossings == @[[15, 14, 5, 4], [3, 2, 13, 12], [9, 8, 11, 10], [7, 6, 1, 0]]
    check result.len == 2
    check result[0].crossings == [[11, 10, 5, 4], [3, 2, 9, 8], [7, 6, 1, 0]]
    check result[0].signs == [-1, -1, -1]
    check result[0].link_components == [newLinkComponent[int](1, 2)]
    check result[1].crossings == [[1, 0, 3, 2]]
    check result[1].signs == [1]
    check result[1].link_components == [newLinkComponent[int](0, 2)]

test "over_or_under_arcs":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check over_or_under_arcs(trefoil, true) == @[@[(0, 2), (1, 3)], @[(1, 2), (2, 3)], @[(2, 2), (0, 3)]]
    check over_or_under_arcs(trefoil, false) == @[@[(0, 3), (1, 2)], @[(1, 3), (2, 2)], @[(2, 3), (0, 2)]]
    # unlinked_hopf has cyclic components lying entirely above/below the rest of the diagram
    let unlinked_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 3, 1, 2]], @[1, -1])
    check over_or_under_arcs(unlinked_hopf, true) == @[@[(0, 2)], @[(1, 2)], @[(0, 1), (1, 3), (0, 1)]]
    check over_or_under_arcs(unlinked_hopf, false) == @[@[(0, 1)], @[(1, 3)], @[(0, 2), (1, 2), (0, 2)]]
    let unoriented_trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]], @[0, 0, 0])
    check over_or_under_arcs(unoriented_trefoil, true) == @[@[(0, 0), (2, 1)], @[(0, 2), (1, 3)], @[(1, 0), (0, 1)]]
    check over_or_under_arcs(unoriented_trefoil, false) == @[@[(0, 1), (2, 0)], @[(0, 3), (1, 2)], @[(1, 1), (0, 0)]]

test "randomize_within_lengths":
    var item = @[@[2], @[1, 3]]
    check randomize_within_lengths(item) == @[@[1, 3], @[2]]

test "pickup_arc_internal":
    var unlinked_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 3, 1, 2]], @[1, -1])
    check pickup_arc_internal(unlinked_hopf, true, @[(0, 1), (1, 3), (0, 1)]) == (@[1, 0], @[])
    # check arc_is_cycle
    # check cross_pos = {1: 0, 0: 1}
    # check num_middle_crossings == 2
    # check link.link_components == @[(0, 1), (0, 2)]
    # check remove_comps == @[true, true]
    # check link.link_components == @[]
    # check elim == @[1, 0]
    # check (crossed_c, crossed_s) == (1, 0)
    # check (start_c, start_s) == (1, 0)
    # check link.crossings == @[[6, 5, 4, 7], [2, 1, 0, 3]]
    # check (incoming_c, incoming_s) == (0, 1)
    # check (outgoing_c, outgoing_s) == (1, 1)
    # check elim_updated == @[1, 0]
    # check link.crossings == @[]
    check unlinked_hopf.crossings == []
    check unlinked_hopf.signs == []
    check unlinked_hopf.unlinked_unknot_components == 2
    check unlinked_hopf.link_components == []

    unlinked_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 5, 1, 2], [3, 4, 4, 5]])
    check pickup_arc_internal(unlinked_hopf, false, @[(1, 2), (0, 2), (1, 2)]) == (@[0, 1], @[])
    # check arc_is_cycle
    # check link.link_components == @[(2, 3), (0, 2)]
    # check remove_comps == @[false, true]
    # check link.link_components == @[(2, 3)]
    check unlinked_hopf.crossings == [[3, 2, 1, 0]]
    check unlinked_hopf.signs == [-1]
    check unlinked_hopf.unlinked_unknot_components == 1
    check unlinked_hopf.link_components == [newLinkComponent[int](0, 3)]

    unlinked_hopf = link_from_PD_code(@[[0, 4, 1, 5], [3, 4, 0, 5], [1, 3, 2, 2]], @[-1, 1, 1])
    check pickup_arc_internal(unlinked_hopf, false, @[(2, 0), (0, 0), (1, 0)]) == (@[0, 1], @[])
    # check not arc_is_cycle
    # check link.link_components == @[(2, 1), (0, 3)]
    # check remove_comps == @[false, true]
    # check link.link_components == @[(2, 1)]
    check unlinked_hopf.crossings == [[1, 0, 3, 2]]
    check unlinked_hopf.signs == [1]
    check unlinked_hopf.unlinked_unknot_components == 1
    check unlinked_hopf.link_components == [newLinkComponent[int](0, 1)]

    let trefoil_4cross = link_from_PD_code(@[[7, 4, 0, 5], [3, 0, 4, 1], [6, 1, 7, 2], [2, 5, 3, 6]])
    check pickup_arc_internal(trefoil_4cross, true, @[(3, 0), (2, 1), (1, 1)]) == (@[2, 1], @[2])
    # check not arc_is_cycle
    # check cross_pos == {2: 0, 1: 1}
    # check num_middle_crossings == 2
    # check elim == @[2, 1]
    # check (crossed_c, crossed_s) == (2, 2)
    # check (start_c, start_s) == (3, 3)
    # check (end_c, end_s) = (0, 0)
    # check (crossed_c, crossed_s) = (1, 2)
    # check (start_c, start_s) == (3, 2)
    # check (end_c, end_s) = (0, 1)
    # check link.crossings == @[[15, 14, 5, 13], [14, 2, 1, 9], [15, 7, 0, 12], [11, 3, 1, 0]]
    # check link.crossings == @[[7, 6, -2, 5], [-2, 3, 1, 0]]
    # check temp_faces == @[@[(0, 0), (1, 2)], @[(0, 1), (1, 1)], @[(0, 3), (1, 3)]]
    # check temp_dual_graph == @[@[(2, (0, 0)), (1, (1, 2))], @[(0, (0, 1)), (2, (1, 1))], @[(1, (0, 3)), (0, (1, 3))]]
    # check distances == @[1, 1, 0]
    # check previous_face == @[2, 2, 0]
    # check last_interface_edge == @[(1, 3), (0, 3), (0, 0)]
    # check shortest_path_edges == @[(0, 3)]
    check trefoil_4cross.crossings == [[7, 6, 9, 8], [11, 10, 1, 0], [3, 2, 5, 4]]
    check trefoil_4cross.signs == [-1, -1, -1]
    check trefoil_4cross.unlinked_unknot_components == 0
    check trefoil_4cross.link_components == [newLinkComponent[int](0, 2)]

test "pickup_arc":
    # neither left or right cycle

    let arc_len_1 = link_from_PD_code(@[[3, 1, 0, 0], [2, 1, 3, 2]])
    check pickup_arc(arc_len_1, true, @[(0, 1)]) == 0
    # check left_cycle_index == -1
    # check right_cycle_index == -1

    # right cycle (with or without left cycle)

    let one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    check pickup_arc(one_vertex_unknot, false, @[(0, 3), (0, 2)]) == 1
    # check left_cycle_index == 1
    # check right_cycle_index == 1
    # check elim == @[]
    # check (cur_c, cur_s) == (0, 2)
    # check new_arc == @[(0, 3)] 
    check one_vertex_unknot.crossings == []
    check one_vertex_unknot.signs == []
    check one_vertex_unknot.unlinked_unknot_components == 1
    check one_vertex_unknot.link_components == []

    let elim_size_2 = link_from_PD_code(@[[5, 4, 0, 5], [3, 0, 4, 1], [6, 1, 7, 2], [7, 3, 6, 2]], @[-1, -1, -1, 1])
    check pickup_arc(elim_size_2, true, @[(0, 2), (1, 3), (2, 3), (3, 1)]) == 3
    # check left_cycle_index == -1
    # check right_cycle_index == 1
    # check elim == @[2, 3]
    # check (cur_c, cur_s) == (1, 3)
    # check new_arc == @[(0, 2)]
    # check new_arc == @[(0, 2)]

    let first_pickup_relabel = link_from_PD_code(@[[5, 4, 0, 5], [6, 1, 7, 2], [3, 0, 4, 1], [7, 3, 6, 2]], @[-1, -1, -1, 1])
    check pickup_arc(first_pickup_relabel, true, @[(0, 2), (2, 3), (1, 3), (3, 1)]) == 3
    # check left_cycle_index == -1
    # check right_cycle_index == 1
    # check elim == @[1, 3]
    # check (cur_c, cur_s) == (1, 3)
    # check new_arc == @[(0, 2)]
    # check new_arc == @[(0, 2)]

    let countdown_len_2 = link_from_PD_code(@[[3, 1, 0, 0], [2, 1, 3, 2]])
    check pickup_arc(countdown_len_2, true, @[(0, 2), (0, 1), (1, 3)]) == 2
    # check left_cycle_index == 1
    # check right_cycle_index == 2
    # check elim == @[]
    # check (cur_c, cur_s) == (1, 3)
    # check new_arc == @[(0, 2), (0, 1)]
    # check new_arc == @[(0, 2), (0, 1)]

    let update_new_arc_1_0 = link_from_PD_code(@[[2, 1, 3, 2], [3, 1, 0, 0]])
    check pickup_arc(update_new_arc_1_0, true, @[(1, 2), (1, 1), (0, 3)]) == 2
    # check left_cycle_index == 1
    # check right_cycle_index == 2
    # check elim == @[]
    # check (cur_c, cur_s) == (0, 3)
    # check new_arc == @[(1, 2), (1, 1)]
    # check new_arc == @[(0, 2), (0, 1)]

    let prev_strand_new_arc_1 = link_from_PD_code(@[[3, 3, 0, 2], [4, 0, 5, 1], [7, 2, 4, 1], [6, 6, 7, 5]])
    check pickup_arc(prev_strand_new_arc_1, true, @[(0, 2), (1, 3), (2, 1), (0, 1)]) == 3
    # check left_cycle_index == 3
    # check right_cycle_index == 3
    # check elim == @[]
    # check (cur_c, cur_s) == (0, 1)
    # check new_arc == @[(0, 2), (1, 3), (2, 1)]
    # check new_arc == @[(2, 1), (1, 3), (2, 1)]

    let update_new_arc_0_0 = link_from_PD_code(@[[3, 0, 4, 1], [5, 4, 0, 5], [6, 1, 7, 2], [7, 3, 6, 2]], @[-1, -1, -1, 1])
    check pickup_arc(update_new_arc_0_0, true, @[(1, 2), (0, 3), (2, 3), (3, 1)]) == 3
    # check left_cycle_index == -1
    # check right_cycle_index == 1
    # check elim == @[2, 3]
    # check (cur_c, cur_s) == (0, 3)
    # check new_arc == @[(1, 2)]
    # check new_arc == @[(0, 2)]

    # left cycle (and no right cycle)

    let left_cycle_no_elim = link_from_PD_code(@[[3, 1, 0, 0], [1, 4, 2, 5], [4, 3, 5, 2]])
    check pickup_arc(left_cycle_no_elim, true, @[(0, 2), (0, 1)]) == 1
    # check left_cycle_index == 1
    # check right_cycle_index == -1
    # check elim == @[]
    # check (cur_c, cur_s) == (0, 1)
    # check new_arc == @[(0, 1)]
    check left_cycle_no_elim.crossings == [[5, 4, 7, 6], [1, 0, 3, 2]]
    check left_cycle_no_elim.signs == [1, 1]
    check left_cycle_no_elim.unlinked_unknot_components == 0
    check left_cycle_no_elim.link_components == [newLinkComponent[int](0, 1), newLinkComponent[int](1, 1)]

    let elim_len_2 = link_from_PD_code(@[[5, 3, 0, 2], [7, 1, 6, 0], [6, 1, 7, 2], [3, 5, 4, 4]], @[1, 1, -1, 1])
    check pickup_arc(elim_len_2, true, @[(0, 2), (1, 1), (2, 3), (0, 1)]) == 3
    # check left_cycle_index == 3
    # check right_cycle_index == -1
    # check elim == @[1, 2]
    # check (cur_c, cur_s) == (0, 1)
    # check new_arc == @[(0, 1)]
    check elim_len_2.crossings == [[1, 0, 3, 2]]
    check elim_len_2.signs == [1]
    check elim_len_2.unlinked_unknot_components == 1
    check elim_len_2.link_components == [newLinkComponent[int](0, 1)]

    let elim_len_2_update = link_from_PD_code(@[[7, 1, 6, 0], [3, 5, 4, 4], [5, 3, 0, 2], [6, 1, 7, 2]])
    check pickup_arc(elim_len_2_update, true, @[(2, 2), (0, 1), (3, 3), (2, 1)]) == 3
    # check left_cycle_index == 3
    # check right_cycle_index == -1
    # check elim == @[0, 3]
    # check (cur_c, cur_s) == (0, 1)
    # check new_arc == @[(0, 1)]
    check elim_len_2_update.crossings == [[1, 0, 3, 2]]
    check elim_len_2_update.signs == [1]
    check elim_len_2_update.unlinked_unknot_components == 1
    check elim_len_2_update.link_components == [newLinkComponent[int](0, 1)]

    let new_arc_len_2 = link_from_PD_code(@[[3, 1, 0, 0], [4, 1, 5, 2], [2, 5, 3, 4]])
    check pickup_arc(new_arc_len_2, true, @[(0, 2), (0, 1), (1, 3)]) == 1
    # check left_cycle_index == 1
    # check right_cycle_index == -1
    # check elim == @[]
    # check (cur_c, cur_s) == (0, 1)
    # check new_arc == @[(0, 1), (1, 3)]
    # check new_arc == @[(0, 2), (1, 3)]
    check new_arc_len_2.crossings == [[7, 6, 5, 4], [3, 2, 1, 0]]
    check new_arc_len_2.signs == [-1, -1]
    check new_arc_len_2.unlinked_unknot_components == 0
    check new_arc_len_2.link_components == [newLinkComponent[int](0, 2), newLinkComponent[int](0, 3)]

    let new_arc_len_2_update = link_from_PD_code(@[[3, 1, 0, 0], [2, 5, 3, 4], [4, 1, 5, 2]])
    check pickup_arc(new_arc_len_2_update, true, @[(0, 2), (0, 1), (2, 3)]) == 1
    # check left_cycle_index == 1
    # check right_cycle_index == -1
    # check elim == @[]
    # check (cur_c, cur_s) == (0, 1)
    # check new_arc == @[(0, 1), (2, 3)]
    # check new_arc == @[(1, 2), (0, 3)]
    check new_arc_len_2_update.crossings == [[7, 6, 5, 4], [3, 2, 1, 0]]
    check new_arc_len_2_update.signs == [-1, -1]
    check new_arc_len_2_update.unlinked_unknot_components == 0
    check new_arc_len_2_update.link_components == [newLinkComponent[int](0, 2), newLinkComponent[int](0, 3)]

test "optimize_overcrossings":
    var L = link_from_PD_code(@[[10, 4, 11, 3], [7, 2, 8, 3], [8, 0, 9, 5],
                                [4, 10, 5, 9], [1, 6, 2, 7], [11, 0, 6, 1]])
    check L.crossings.len == 6
    check L.simplify(SimplificationMode.level) == 0
    check L.optimize_overcrossings() == 1

test "untwist_diagram_once":
    let one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    check untwist_diagram_once(one_vertex_unknot) == (@[0], @[])
    check one_vertex_unknot.crossings == []
    check one_vertex_unknot.signs == []
    check one_vertex_unknot.unlinked_unknot_components == 1
    check one_vertex_unknot.link_components == []

    let right_degenerate = link_from_PD_code(@[[3, 4, 0, 5], [4, 1, 5, 0], [1, 2, 2, 3]])
    check untwist_diagram_once(right_degenerate) == (@[2], @[0, 1])
    check right_degenerate.crossings == [[5, 4, 7, 6], [1, 0, 3, 2]]
    check right_degenerate.signs == [1, 1]
    check right_degenerate.unlinked_unknot_components == 0
    check right_degenerate.link_components == [newLinkComponent[int](0, 1), newLinkComponent[int](0, 2)]

    let left_degenerate = link_from_PD_code(@[[3, 2, 0, 3], [5, 0, 4, 1], [1, 4, 2, 5]])
    check untwist_diagram_once(left_degenerate) == (@[0], @[0, 1])
    check left_degenerate.crossings == [[7, 6, 5, 4], [3, 2, 1, 0]]
    check left_degenerate.signs == [-1, -1]
    check left_degenerate.unlinked_unknot_components == 0
    check left_degenerate.link_components == [newLinkComponent[int](0, 2), newLinkComponent[int](1, 2)]

    # TODO: nondegenerate_no_move_comp

    let nondegenerate_move_comp = link_from_PD_code(@[[5, 2, 0, 3], [1, 0, 2, 1], [3, 4, 4, 5]])
    check untwist_diagram_once(nondegenerate_move_comp) == (@[0], @[0, 1])
    check nondegenerate_move_comp.crossings == [[6, 2, 1, 5], [7, 3, 0, 4]]
    check nondegenerate_move_comp.signs == [-1, -1]
    check nondegenerate_move_comp.unlinked_unknot_components == 0
    check nondegenerate_move_comp.link_components == [newLinkComponent[int](0, 3)]

    let middle_s_1 = link_from_PD_code(@[[5, 3, 0, 2], [1, 1, 2, 0], [3, 5, 4, 4]])
    check untwist_diagram_once(middle_s_1) == (@[0], @[0, 1])
    check middle_s_1.crossings == [[6, 7, 3, 2], [5, 4, 0, 1]]
    check middle_s_1.signs == [1, 1]
    check middle_s_1.unlinked_unknot_components == 0
    check middle_s_1.link_components == [newLinkComponent[int](1, 2)]

    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check untwist_diagram_once(trefoil) == (@[], @[])

test "untwist_diagram":
    # a connect sum of opposite-handed trefoils with three extra twists
    let L = link_from_PD_code(@[[2, 0, 3, 17], [16, 2, 17, 1], [0, 16, 1, 15],
                                [9, 6, 10, 7], [7, 10, 8, 11], [11, 8, 12, 9],
                                [3, 15, 4, 14], [13, 5, 14, 4], [5, 13, 6, 12]])
    check untwist_diagram(L) == 3

test "pickup_simplify":
    # 13n3370, with one crossing change via band attachment, with another band attached (and thus the diagram of an unknot)
    var cc2 = link_from_PD_code(@[[0, 21, 1, 22], [20, 1, 21, 2], [2, 19, 36, 20], [18, 32, 19, 13], [17, 30, 18, 27], [25, 16, 26, 17], [4, 16, 5, 15], [14, 4, 15, 3], [7, 13, 8, 12], [11, 28, 12, 29], [10, 24, 11, 23], [22, 10, 23, 9], [29, 8, 0, 9], [6, 27, 7, 28], [24, 5, 25, 6], [26, 14, 33, 31], [30, 31, 34, 32], [33, 3, 37, 35], [35, 37, 36, 34]])
    check pickup_simplify(cc2, 1) == 19

    # check that twists are undone
    let three_rings_with_twist = link_from_PD_code(@[[5, 7, 0, 6], [6, 0, 7, 1], [1, 5, 2, 4], [9, 2, 8, 3], [3, 8, 4, 9]])
    check pickup_simplify(three_rings_with_twist, 0) == 1

test "reverse_type_i":
    var one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    reverse_type_i(one_vertex_unknot, (0, 2), false)
    # check given_strand_sign == 1
    # check given_strand_index == 1
    # check new_crossing_sign == -1
    check one_vertex_unknot.crossings == [[3, 6, 5, 0], [7, 2, 1, 4]]
    check one_vertex_unknot.signs == [-1, -1]

    one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]], @[0])
    reverse_type_i(one_vertex_unknot, (0, 2), false)
    # check given_strand_sign == 0
    # check given_strand_index == 1
    # check new_crossing_sign == 0
    check one_vertex_unknot.crossings == [[3, 6, 5, 0], [7, 2, 1, 4]]
    check one_vertex_unknot.signs == [0, 0]

    one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    reverse_type_i(one_vertex_unknot, (0, 0), false)
    # check given_strand_sign == -1
    # check given_strand_index == 3
    # check new_crossing_sign == -1
    check one_vertex_unknot.crossings == [[7, 2, 1, 4], [3, 6, 5, 0]]
    check one_vertex_unknot.signs == [-1, -1]

    one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    reverse_type_i(one_vertex_unknot, (0, 2), true)
    # check given_strand_sign == 1
    # check given_strand_index == 0
    # check new_crossing_sign == 1
    check one_vertex_unknot.crossings == [[3, 5, 4, 0], [2, 1, 7, 6]]
    check one_vertex_unknot.signs == [-1, 1]

    one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]], @[0])
    reverse_type_i(one_vertex_unknot, (0, 2), true)
    # check given_strand_sign == 0
    # check given_strand_index == 0
    # check new_crossing_sign == 0
    check one_vertex_unknot.crossings == [[3, 5, 4, 0], [2, 1, 7, 6]]
    check one_vertex_unknot.signs == [0, 0]

    one_vertex_unknot = link_from_PD_code(@[[0, 1, 1, 0]])
    reverse_type_i(one_vertex_unknot, (0, 0), true)
    # check given_strand_sign == -1
    # check given_strand_index == 2
    # check new_crossing_sign == 1
    check one_vertex_unknot.crossings == [[6, 2, 1, 7], [5, 4, 0, 3]]
    check one_vertex_unknot.signs == [-1, 1]

test "reverse_type_ii":
    var same_orientation = link_from_PD_code(@[[3, 3, 0, 2], [1, 1, 2, 0]])
    reverse_type_ii(same_orientation, (0, 2), (1, 2))
    check same_orientation.crossings == [[1, 0, 8, 11], [5, 4, 15, 14], [2, 13, 12, 3], [10, 9, 7, 6]]
    check same_orientation.signs == [1, 1, -1, 1]

    same_orientation = link_from_PD_code(@[[3, 3, 0, 2], [1, 1, 2, 0]])
    reverse_type_ii(same_orientation, (1, 2), (0, 2))
    check same_orientation.crossings == [[1, 0, 15, 14], [5, 4, 8, 11], [6, 13, 12, 7], [10, 9, 3, 2]]
    check same_orientation.signs == [1, 1, -1, 1]
