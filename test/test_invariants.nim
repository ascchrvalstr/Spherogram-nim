import std/unittest
import std/json
import std/tables
import std/strscans

import ../src/links
import ../src/invariants

test "linking_matrix":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.linking_matrix() == @[@[0]]
    let l2a1_0 = link_from_PD_code(@[[4, 1, 3, 2], [2, 3, 1, 4]])
    check l2a1_0.linking_matrix() == @[@[0, -1], @[-1, 0]]
    let l2a1_1 = link_from_PD_code(@[[4, 2, 3, 1], [2, 4, 1, 3]])
    check l2a1_1.linking_matrix() == @[@[0, 1], @[1, 0]]

test "knot_floer_homology":
    const hfk_data_json = staticRead"HFK_data.json"
    let knot_list = parseJson(hfk_data_json)
    for knot_node in knot_list.items():
        var pd_code = newSeq[Crossing]()
        for pd_code_item in knot_node["PD_code"].items():
            assert pd_code_item.len == 4
            var pd_code_crossing: array[0..3, int]
            for s in 0 ..< 4:
                pd_code_crossing[s] = pd_code_item[s].getInt()
            pd_code.add(pd_code_crossing)
        var hfk = link_from_PD_code(pd_code).knot_floer_homology()
        check hfk.modulus == knot_node["modulus"].getInt()
        var ranks = knot_node["ranks"]
        check hfk.ranks.len == ranks.len
        # echo hfk.ranks.len, " ", ranks.len
        for (key, val) in ranks.pairs():
            var a, b: int
            discard scanf(key, "($i, $i)", a, b)
            check hfk.ranks[(a, b)] == val.getInt()
            # echo a, " ", b, " ", val.getInt()
        check hfk.total_rank == knot_node["total_rank"].getInt()
        check hfk.seifert_genus == knot_node["seifert_genus"].getInt()
        check hfk.fibered == knot_node["fibered"].getBool()
        check hfk.L_space_knot == knot_node["L_space_knot"].getBool()
        check hfk.tau == knot_node["tau"].getInt()
        check hfk.nu == knot_node["nu"].getInt()
        check hfk.epsilon == knot_node["epsilon"].getInt()

test "white_graph":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.white_graph() == @[@[(1, (0, -1)), (1, (1, -1)), (1, (2, -1))], @[(0, (0, -1)), (0, (1, -1)), (0, (2, -1))]]
    # check strand_to_face == @[[0, 1, 2, 3], [2, 1, 4, 3], [4, 1, 0, 3]]
    # check edges == @[@[(2, (0, 1)), (4, (2, 1))], @[(3, (0, -1)), (3, (1, -1)), (3, (2, -1))], @[(0, (0, 1)), (4, (1, 1))], @[(1, (0, -1)), (1, (1, -1)), (1, (2, -1))], @[(2, (1, 1)), (0, (2, 1))]]
    # check conn_comps == @[@[0, 2, 4], @[1, 3]]
    # check conn_comps == @[@[0, 2, 4], @[1, 3]]
    # check index_within_comp_1 == @[-1, 0, -1, 1, -1]

test "goeritz_matrix":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.goeritz_matrix()[0] == @[@[3]]
    # check m == @[@[0, -3], @[-3, 0]]
    # check m == @[@[3, -3], @[-3, 3]]

test "signature":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.signature() == 2

    let figure_eight = link_from_PD_code(@[[7, 4, 0, 5], [3, 0, 4, 1], [1, 7, 2, 6], [5, 3, 6, 2]])
    check figure_eight.signature() == 0

    let knot_51 = link_from_PD_code(@[[9, 4, 0, 5], [5, 0, 6, 1], [1, 6, 2, 7], [7, 2, 8, 3], [3, 8, 4, 9]])
    check knot_51.signature() == 4

    let knot_52 = link_from_PD_code(@[[4, 0, 5, 9], [0, 6, 1, 5], [8, 2, 9, 1], [2, 8, 3, 7], [6, 4, 7, 3]])
    check knot_52.signature() == -2

test "colorability_matrix":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.colorability_matrix() == @[@[-1, -1, 2], @[2, -1, -1], @[-1, 2, -1]]
    let hopf_positive = link_from_PD_code(@[[2, 1, 3, 0], [1, 2, 0, 3]])
    check hopf_positive.colorability_matrix() == @[@[-2, 2], @[2, -2]]
    let hopf_negative = link_from_PD_code(@[[0, 2, 1, 3], [3, 1, 2, 0]])
    check hopf_negative.colorability_matrix() == @[@[-2, 2], @[2, -2]]

test "determinant":
    let trefoil = link_from_PD_code(@[[5, 2, 0, 3], [3, 0, 4, 1], [1, 4, 2, 5]])
    check trefoil.determinant("color") == 3
    check trefoil.determinant("goeritz") == 3
    let hopf_positive = link_from_PD_code(@[[2, 1, 3, 0], [1, 2, 0, 3]])
    check hopf_positive.determinant("color") == 2
    check hopf_positive.determinant("goeritz") == 2
    let figure_eight = link_from_PD_code(@[[7, 4, 0, 5], [3, 0, 4, 1], [1, 7, 2, 6], [5, 3, 6, 2]])
    check figure_eight.determinant("color") == 5
    check figure_eight.determinant("goeritz") == 5
    let knot_51 = link_from_PD_code(@[[9, 4, 0, 5], [5, 0, 6, 1], [1, 6, 2, 7], [7, 2, 8, 3], [3, 8, 4, 9]])
    check knot_51.determinant("color") == 5
    check knot_51.determinant("goeritz") == 5
