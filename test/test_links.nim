import std/unittest
import std/sequtils

import ../src/links

# performance:
# construction of knot 13_1 from PD code: ~100000 per second vs ~3500 for Spherogram

test "writhe":
    let left_handed_trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check left_handed_trefoil.writhe() == -3

test "rotate_crossing":
    var left_handed_trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    left_handed_trefoil.rotate_crossing(0, 0)
    check left_handed_trefoil.crossings == [[11, 10, 5, 4], [3, 2, 9, 8], [7, 6, 1, 0]]
    check left_handed_trefoil.signs == [-1, -1, -1]
    check left_handed_trefoil.link_components == [newLinkComponent[int](0, 2)]

    left_handed_trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    left_handed_trefoil.rotate_crossing(0, 3)
    check left_handed_trefoil.crossings == [[10, 5, 4, 11], [2, 1, 9, 8], [7, 6, 0, 3]]
    check left_handed_trefoil.signs == [1, -1, -1]
    check left_handed_trefoil.link_components == [newLinkComponent[int](0, 1)]

    left_handed_trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]], @[0, 0, 0])
    left_handed_trefoil.rotate_crossing(0, 2)
    check left_handed_trefoil.crossings == [[5, 4, 11, 10], [1, 0, 9, 8], [7, 6, 3, 2]]
    check left_handed_trefoil.signs == [0, 0, 0]
    check left_handed_trefoil.link_components == [newLinkComponent[int](0, 2)]

    var two_vertex_unknot = link_from_PD_code(@[[3, 1, 0, 0], [1, 3, 2, 2]])
    two_vertex_unknot.rotate_crossing(0, 1)
    check two_vertex_unknot.crossings == [[3, 5, 4, 0], [2, 1, 7, 6]]
    check two_vertex_unknot.signs == [-1, 1]
    check two_vertex_unknot.link_components == [newLinkComponent[int](0, 2)]

test "is_alternating":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.is_alternating()
    let silly_hopf_link = link_from_PD_code(@[[4, 8, 1, 5], [3, 6, 4, 5], [6, 3, 7, 2], [1, 8, 2, 7]])
    check not silly_hopf_link.is_alternating()

test "split_link_diagram":
    let trefoil_and_unknot = link_from_PD_code(@[[5, 2, 0, 3], [6, 8, 7, 7], [3, 0, 4, 1], [9, 8, 6, 9], [1, 4, 2, 5]])
    let conn_comps = trefoil_and_unknot.split_link_diagram()
    check conn_comps[0].crossings == [[11, 10, 5, 4], [3, 2, 9, 8], [7, 6, 1, 0]]
    check conn_comps[0].signs == [-1, -1, -1]
    check conn_comps[0].unlinked_unknot_components == 0
    check conn_comps[0].link_components == [newLinkComponent[int](0, 2)]
    check conn_comps[1].crossings == [[6, 5, 3, 2], [7, 1, 0, 4]]
    check conn_comps[1].signs == [1, -1]
    check conn_comps[1].unlinked_unknot_components == 0
    check conn_comps[1].link_components == [newLinkComponent[int](0, 1)]

test "disjoint_union":
    var knot61 = link_from_PD_code(@[[6, 11, 7, 0], [0, 5, 1, 6], [10, 2, 11, 1], [2, 10, 3, 9], [8, 4, 9, 3], [4, 8, 5, 7]], unlinked_unknot_components = 1)
    var knot41 = link_from_PD_code(@[[7, 4, 0, 5], [3, 0, 4, 1], [1, 7, 2, 6], [5, 3, 6, 2]], unlinked_unknot_components = 2)
    var union_61_41 = disjoint_union(knot61, knot41)
    check (union_61_41.PD_code() ==
           @[[11, 4, 0, 5],     # 6_1
             [5, 10, 6, 11],
             [3, 7, 4, 6],
             [7, 3, 8, 2],
             [1, 9, 2, 8],
             [9, 1, 10, 0],
             [19, 16, 12, 17],  # 4_1 with
             [15, 12, 16, 13],  # shifted
             [13, 19, 14, 18],  # edge labels
             [17, 15, 18, 14]])
    check (union_61_41.signs ==
           @[-1, -1, 1, 1, 1, 1,  # 6_1
             -1, -1, 1, 1])       # 4_1
    # test that the numbers of unlinked unknot
    # components have been correctly added
    check union_61_41.unlinked_unknot_components == 3
    check union_61_41.link_components == @[newLinkComponent[int](0, 2),
                                           newLinkComponent[int](6, 2)]

test "connected_sum":
    let link1 = link_from_PD_code(@[[0, 1, 1, 2], [3, 3, 0, 2]])
    let link2 = link_from_PD_code(@[[0, 0, 1, 1]])
    var conn_sum = connected_sum(link1, link2)
    check conn_sum.crossings == [[9, 2, 1, 7], [5, 4, 8, 3], [6, 0, 11, 10]]
    check conn_sum.signs == [-1, 1, 1]
    check conn_sum.link_components == [newLinkComponent[int](0, 2)]
    let link3 = link_from_PD_code(@[[0, 0, 1, 1], [2, 2, 3, 3]])
    conn_sum = connected_sum(link1, link3)
    check conn_sum.link_components == [newLinkComponent[int](0, 2), newLinkComponent[int](3, 1)]

test "sublink":
    let trefoil_linked_unknot = link_from_PD_code(@[[9, 7, 8, 0], [0, 8, 1, 9], [4, 1, 5, 2], [2, 5, 3, 6], [6, 3, 7, 4]])
    let unknot = trefoil_linked_unknot.sublink(@[true, false])
    # check strand_to_comp == @[[0, -1, 0, -1], [-1, 0, -1, 0], [-1, -1, -1, -1], [-1, -1, -1, -1], [-1, -1, -1, -1]]
    # check new_crossing_labels == @[-1, -1, -1, -1, -1]
    # check unlinked == @[true, false]
    check unknot.crossings == []
    check unknot.signs == []
    check unknot.unlinked_unknot_components == 1
    check unknot.link_components == []
    let trefoil = trefoil_linked_unknot.sublink(@[false, true])
    # check strand_to_comp == @[[-1, 1, -1, 1], [1, -1, 1, -1], [1, 1, 1, 1], [1, 1, 1, 1], [1, 1, 1, 1]]
    # check new_crossing_labels == @[-1, -1, 0, 1, 2]
    # check unlinked == @[false, false]
    check trefoil.crossings == @[[11, 10, 5, 4], [3, 2, 9, 8], [7, 6, 1, 0]]
    check trefoil.signs == @[-1, -1, -1]
    check trefoil.unlinked_unknot_components == 0
    check trefoil.link_components == @[newLinkComponent[int](2, 2)]

test "linking_number":
    let unlinked_hopf = link_from_PD_code(@[[2, 1, 3, 0], [3, 1, 2, 0]])
    check unlinked_hopf.linking_number() == 0
    let hopf_positive = link_from_PD_code(@[[2, 1, 3, 0], [1, 2, 0, 3]])
    check hopf_positive.linking_number() == 1
    let hopf_negative = link_from_PD_code(@[[0, 2, 1, 3], [3, 1, 2, 0]])
    check hopf_negative.linking_number() == -1
    # check that we are not counting crossings within the same link component
    let trefoil_linked_unknot = link_from_PD_code(@[[9, 7, 8, 0], [0, 8, 1, 9], [4, 1, 5, 2], [2, 5, 3, 6], [6, 3, 7, 4]])
    check trefoil_linked_unknot.linking_number() == -1
    # performance:
    # trefoil_linked_unknot.linking_number(): ~5000000 per second vs ~350000 for Spherogram

test "is_planar":
    let bad = link_from_PD_code(@[[0, 1, 0, 1]])
    check not bad.is_planar()
    let L = link_from_PD_code(@[[1, 7, 2, 6], [7, 4, 8, 5], [3, 8, 0, 9], [5, 3, 6, 2], [9, 0, 4, 1]])
    check L.is_planar()
    let S = link_from_PD_code(@[[1, 1, 2, 2], [3, 3, 4, 4]])
    check S.is_planar()
    let N = link_from_PD_code(@[[0, 0, 1, 1], [2, 3, 2, 3]])
    check not N.is_planar()
    check mapIt(N.split_link_diagram(), it.is_planar()) == [true, false]

test "pieces":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.pieces() == @[@[(0, 2), (1, 1), (1, 3), (2, 0)], @[(1, 2), (2, 1), (2, 3), (0, 0)], @[(2, 2), (0, 1), (0, 3), (1, 0)]]
    let unoriented_trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]], @[0, 0, 0])
    check unoriented_trefoil.pieces() == @[@[(0, 0), (2, 3), (2, 1), (1, 2)], @[(0, 2), (1, 1), (1, 3), (2, 0)], @[(1, 0), (0, 3), (0, 1), (2, 2)]]
    let unlinked_hopf = link_from_PD_code(@[[1, 3, 0, 2], [0, 3, 1, 2]], @[1, -1])
    check unlinked_hopf.pieces() == @[@[(0, 2), (1, 0)], @[(1, 2), (0, 0)]]
