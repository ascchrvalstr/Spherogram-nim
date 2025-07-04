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
