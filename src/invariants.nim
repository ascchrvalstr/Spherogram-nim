import std/sequtils
import std/enumerate
import std/deques
import std/strformat
import std/algorithm

import manu

import links
import knot_floer_homology

proc linking_matrix*[T](link: Link[T]): seq[seq[int]] =
    let num_crossings = link.crossings.len
    var strand_to_comp = newSeqWith(num_crossings, [-1, -1, -1, -1])
    for comp_index, comp in enumerate(link.link_components):
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        while strand_to_comp[cur_c][cur_s] == -1:
            strand_to_comp[cur_c][cur_s] = comp_index
            strand_to_comp[cur_c][(cur_s+2) mod 4] = comp_index
            (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
    let num_components = link.link_components.len
    var linking_matrix = newSeqWith(num_components, newSeqWith(num_components, 0))
    for (strands, sign) in zip(strand_to_comp, link.signs):
        if strands[0] != strands[1]:
            linking_matrix[strands[0]][strands[1]] += sign
            linking_matrix[strands[1]][strands[0]] += sign
    for i in 0 ..< num_components:
        for j in 0 ..< num_components:
            linking_matrix[i][j] = linking_matrix[i][j] div 2
    return linking_matrix

proc knot_floer_homology*[T](link: Link[T], prime: int = 2): HFK =
    if link.link_components.len + link.unlinked_unknot_components > 1:
        raise newException(ValueError, &"knot Floer homology only works for knots, but this has {link.link_components.len + link.unlinked_unknot_components} components")
    if link.link_components.len == 0 and link.unlinked_unknot_components == 1:
        return link_from_PD_code(@[[4, 2, 5, 1], [2, 0, 3, 5], [3, 0, 4, 1]]).knot_floer_homology()
    return pd_code_to_hfk(repr(link.PD_code()), prime)

proc add_edge*[T](edges: var seq[seq[(int, T)]], u: int, v: int, edge_data: T): void =
    edges[u].add((v, edge_data))
    if u != v:
        edges[v].add((u, edge_data))

proc white_graph*[T](link: Link[T]): seq[seq[(int, (int, int))]] =
    discard """
    Return the white graph of a non-split link projection.

    This method generates a multigraph whose vertices correspond
    to the faces of the diagram, with an edge joining two
    vertices whenever the corresponding faces contain opposite
    corners at some crossing. The vertex corresponding to a face
    is the index of the face in the list returned by Link.faces().

    According to the conventions of "Gordon, C. McA. and
    Litherland, R. A, 'On the signature of a link', Inventiones
    math. 47, 23-69 (1978)", in a checkerboard coloring of a link
    diagram the unbounded region is always the first white region.
    Of course, the choice of which region is unbounded is
    arbitrary; it is just a matter of which region on S^2 contains
    the point at infinity.  In this method an equivalent arbitrary
    choice is made by just returning the second component of the
    multigraph, as determined by a homegrown re-implementation of
    Graph.connected_components(). (By the documentation of
    Graph.connected_components(), the second component is no
    bigger than the first.)

    Note that this may produce a meaningless result in the case of
    a split link diagram.  Consequently if the diagram is split,
    i.e if the multigraph has more than 2 components, a ValueError
    is raised::

        nim: K = Link('5_1')
        nim: K.white_graph()
        Subgraph of (): Multi-graph on 2 vertices

    WARNING: While there is also a "black_graph" method, it need
    not be the case that these two graphs are complementary in the
    expected way.
    """
    let num_crossings = link.crossings.len
    var strand_to_face = newSeq[array[0..3, int]](num_crossings)
    let faces = link.faces()
    let num_faces = faces.len
    for face_index, face in enumerate(faces):
        for strand in face:
            strand_to_face[strand.crossing][strand.strand_index] = face_index
    # echo strand_to_face
    var edges = newSeq[seq[(int, (int, int))]](num_faces)
    for c in 0 ..< num_crossings:
        edges.add_edge(strand_to_face[c][0], strand_to_face[c][2], (c, 1))
        edges.add_edge(strand_to_face[c][1], strand_to_face[c][3], (c, -1))
    # echo edges
    var conn_comps = newSeq[seq[int]]()
    var visited = newSeq[bool](num_faces)
    for face in 0 ..< num_faces:
        if not visited[face]:
            var bfs_queue = toDeque([face])
            visited[face] = true
            var cur_conn_comp = @[face]
            while bfs_queue.len > 0:
                var cur_face = bfs_queue.popFirst()
                for (next_face, _) in edges[cur_face]:
                    if not visited[next_face]:
                        visited[next_face] = true
                        cur_conn_comp.add(next_face)
                        bfs_queue.addLast(next_face)
            sort(cur_conn_comp)
            conn_comps.add(cur_conn_comp)
    # echo conn_comps
    if conn_comps.len > 2:
        raise newException(ValueError, &"the link diagram is split because the pre-white graph has {conn_comps.len} > 2 connected components")
    sort(conn_comps, proc (x, y: seq[int]): int = y.len - x.len)
    # echo conn_comps
    var index_within_comp_1 = newSeqWith(num_faces, -1)
    for (face_index, face) in enumerate(conn_comps[1]):
        index_within_comp_1[face] = face_index
    # echo index_within_comp_1
    var new_graph = newSeq[seq[(int, (int, int))]](conn_comps[1].len)
    for face_u in conn_comps[1]:
        for (face_v, edge_data) in edges[face_u]:
            if index_within_comp_1[face_v] != -1:
                new_graph[index_within_comp_1[face_u]].add((index_within_comp_1[face_v], edge_data))
    return new_graph

proc goeritz_matrix*[T](link: Link[T]): (seq[seq[int]], seq[seq[(int, (int, int))]]) =
    var G = white_graph(link)
    let N = G.len
    var m = newSeqWith(N, newSeq[int](N))
    for i in 0 ..< N:
        for (j, edge_data) in G[i]:
            if i <= j:
                m[i][j] += edge_data[1]
                m[j][i] = m[i][j]
    # echo m
    for j in 0 ..< N:
        var column_sum = 0
        for i in 0 ..< N:
            column_sum += m[i][j]
        m[j][j] = -column_sum
    # echo m
    var m_remove_r0_c0 = newSeqWith(N-1, newSeq[int](N-1))
    for i in 0 ..< N-1:
        for j in 0 ..< N-1:
            m_remove_r0_c0[i][j] = m[i+1][j+1]
    return (m_remove_r0_c0, G)

proc signature*[T](link: Link[T]): int =
    if link.crossings.len == 0 and link.unlinked_unknot_components == 1:
        return 0
    let (m0, G) = goeritz_matrix(link)
    let size = m0.len
    var m = newSeqWith(size, newSeq[float64](size))
    for i in 0 ..< size:
        for j in 0 ..< size:
            m[i][j] = float64(m0[i][j])
    var ans = 0
    for eigenvalue in matrix(m).eig().getRealEigenvalues():
        if eigenvalue > 0:
            ans += 1
        elif eigenvalue < 0:
            ans -= 1
    for i in 0 ..< G.len:
        for (j, edge_data) in G[i]:
            if i <= j and edge_data[1] == link.signs[edge_data[0]]:
                ans += edge_data[1]
    return -ans

proc colorability_matrix*[T](link: Link[T]): seq[seq[int]] =
    let edges = link.pieces()
    var m = newSeqWith(link.crossings.len, newSeq[int](edges.len))
    for (ind, s) in enumerate(edges):
        for (c, i) in s:
            if i mod 2 == 1:
                m[c][ind] += 1
            else:
                m[c][ind] -= 1
    return m

proc determinant*[T](link: Link[T], method0: string = "goeritz"): int =
    if link.crossings.len == 0 and link.unlinked_unknot_components == 1:
        return 1
    # "method" is a reserved keyword in Nim, so the parameter has to be renamed to "method0"
    if method0 == "color":
        let M = link.colorability_matrix()
        let size = link.crossings.len-1
        var N = newSeqWith(size, newSeq[float64](size))
        for i in 0 ..< size:
            for j in 0 ..< size:
                N[i][j] = float64(M[i+1][j+1])
        # echo matrix(N).det()
        # note that float to int conversion in Nim is similar to C++ in that the number is rounded towards zero, so we need to add 0.5 and then convert to int
        return int(abs(matrix(N).det()) + 0.5)
    if method0 == "goeritz":
        let goeritz_matrix = link.goeritz_matrix()[0]
        let size = goeritz_matrix.len
        var N = newSeqWith(size, newSeq[float64](size))
        for i in 0 ..< size:
            for j in 0 ..< size:
                N[i][j] = float64(goeritz_matrix[i][j])
        # echo matrix(N).det()
        return int(abs(matrix(N).det()) + 0.5)
    raise newException(AssertionDefect, "Not implemented")
