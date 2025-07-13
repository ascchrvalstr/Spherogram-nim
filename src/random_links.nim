discard """
Generating random knots starting with randomly generated planar maps.

Code written as part of a 2013 Illinois Geometry Lab project by

Nathan Dunfield, Alexander Ehrenberg, Shiladitya Bhattacharyya, and Dongming Lei

Details to hopefully appear in some paper or other. I wouldn't hold
your breath, though.
"""
import std/options
import std/strformat
import std/random
import std/sequtils
import std/enumerate
import std/algorithm
import std/sugar

import planarmap
import links
import simplify
import twist

type LinkGenerationError* = object of CatchableError

proc random_link_internal*[T](num_crossings: int, edge_conn_param: int = 4, num_link_comps: int = 0, max_tries: int = 100): Link[T] =
    discard """
    Returns a dictionary of endpoints of edges in the form:

    (signed edge) -> (vertex, position)
    """
    if num_crossings < 0:
        raise newException(ValueError, &"the number of crossings {num_crossings} cannot be negative")
    if num_link_comps < 0:
        raise newException(ValueError, &"the number of link components {num_link_comps} cannot be negative")
    if max_tries <= 0:
        raise newException(ValueError, &"the maximum number of attempts {max_tries} cannot be negative or zero")
    if num_crossings == 0:
        if num_link_comps > 0:
            raise newException(ValueError, &"a zero-crossing link diagram is empty, so it cannot admit {num_link_comps} link components")
        return link_from_PD_code_with_extra_info[T](@[])
    if num_crossings == 1:
        if num_link_comps > 1:
            raise newException(ValueError, &"a one-crossing link diagram must be an unknot after a reverse Reidemeister I move, so it must have one link component, so it cannot have {num_link_comps} link components")
        # there are exactly two possibilities depending on whether the opposite of strand 0 is strand 1 or strand 3
        if rand(0..1) == 0:
            return link_from_PD_code_with_extra_info[T](@[[0, 0, 1, 1]])
        else:
            return link_from_PD_code_with_extra_info[T](@[[0, 1, 1, 0]])
    var data0 = raw_random_map(num_crossings, edge_conn_param, num_link_comps, max_tries)
    if data0.isNone:
        raise newException(LinkGenerationError, &"did not generate the requested link after {max_tries} tries")
    var data = data0.get()
    var vertex_adjacencies = newSeq[seq[int]]()
    for (vertex, adjacencies) in data:
        if rand(0..1) == 1:
            vertex_adjacencies.add(concat(adjacencies[1 ..< adjacencies.len], adjacencies[0 ..< 1]))
        else:
            vertex_adjacencies.add(adjacencies)
    # echo vertex_adjacencies
    var edges_u = newSeq[(int, int)](num_crossings*2)
    var edges_v = newSeq[(int, int)](num_crossings*2)
    for (v, edges) in enumerate(vertex_adjacencies):
        for (i, edge) in enumerate(edges):
            if edge > 0:
                edges_u[edge-1] = (v, i)
            else:
                edges_v[-edge-1] = (v, i)
    # echo edges_u
    # echo edges_v
    var pd_code = newSeq[Crossing](num_crossings)
    for (k, strand_pair) in enumerate(zip(edges_u, edges_v)):
        let (strand1, strand2) = strand_pair
        let (a, i) = strand1
        let (b, j) = strand2
        pd_code[a][i] = k
        pd_code[b][j] = k
    return link_from_PD_code_with_extra_info[T](pd_code)

proc longest_components*[T](link: Link[T], num_requested_components: int): seq[int] =
    let num_components = link.link_components.len
    if num_requested_components < 0 or num_requested_components > num_components:
        raise newException(ValueError, &"num_requested_components = {num_requested_components} is not within [0, num_components] = [0, {num_components}]")
    let num_crossings = link.crossings.len
    var strand_to_comp = newSeqWith(num_crossings, [-1, -1, -1, -1])
    for comp_index, comp in enumerate(link.link_components):
        var cur_c = comp.crossing
        var cur_s = comp.strand_index
        while strand_to_comp[cur_c][cur_s] == -1:
            strand_to_comp[cur_c][cur_s] = comp_index
            strand_to_comp[cur_c][(cur_s+2) mod 4] = comp_index
            (cur_c, cur_s) = link.next_strand((cur_c, cur_s))
    var components_by_self_crossings = newSeq[(int, int)](num_components)
    for i in 0 ..< num_components:
        components_by_self_crossings[i][1] = i
    for c in 0 ..< num_crossings:
        if strand_to_comp[c][0] == strand_to_comp[c][1]:
            components_by_self_crossings[strand_to_comp[c][0]][0] += 1
    sort(components_by_self_crossings, Descending)
    # echo components_by_self_crossings
    return collect:
        for i in 0 ..< num_requested_components:
            components_by_self_crossings[i][1]

proc random_link_with_extra_info*[T](
    crossings: int, num_components0: Option[int] = none(int),
    initial_map_gives_link: bool = false, alternating: bool = false,
    consistent_twist_regions: bool = false,
    simplification_mode: Option[SimplificationMode] = some(SimplificationMode.basic),
    prime_decomposition: bool = true, return_all_pieces: bool = false,
    max_tries: int = 100
): seq[Link[T]] =
    discard """
    Generates a random link from a model that starts with a random
    4-valent planar graph sampled with the uniform distribution by
    Schaeffer's `PlanarMap program.
    <http://www.lix.polytechnique.fr/~schaeffe/PagesWeb/PlanarMap/index-en.html>`_

    The ``crossings`` argument specifies the number of vertices of the
    initial planar graph G; the number of crossing in the returned knot
    will typically be less. The meanings of the optional arguments are as
    follows:

    1. ``num_components``: The number of components of the returned link.
       The link naively associated to G may have too few or too many
       components. The former situation is resolved by picking another G,
       and the latter by either

       a. Taking the sublink consisting of the components with the largest
          self-crossing numbers.

       b. Resampling G until the desired number of components is achieved;
          this can take a very long time as the expected number of
          components associated to G grows linearly in the number of
          vertices.

       When the argument ``initial_map_gives_link`` is ``False`` the
       program does (a) and this is the default behavior. If you want (b)
       set this argument to ``True``.

       To get the entire link associated to G, set ``num_components`` to
       ```any```, which is also the default.

    2. The 4-valent vertices of G are turned into crossings by flipping a
       fair coin. If you want the unique alternating diagram associated to
       G, pass ``alternating=True``.  If you want there to be no
       obvious Type II Reidemeister moves, pass
       ``consistent_twist_regions=True``.

    3. ``simplify``: Whether and how to try to reduce the number of
       crossings of the link via Reidemeister moves using the method
       ``Link.simplify``.  For no simplification, set ``simplify=None``;
       otherwise set ``simplify`` to be the appropriate mode for
       ``Link.simplify``, for example ``basic`` (the default), ``level``,
       or ``global``.

    4. ``prime_decomposition``:  The initial link generated from G may not
       be prime (and typically isn't if ``initial_map_gives_link`` is
       ``False``). When set (the default), the program undoes any connect
       sums that are "diagrammatic obvious", simplifies the result, and
       repeats until pieces are "diagrammatically prime".  If
       ``return_all_pieces`` is ``False`` (the default) then only the
       largest (apparently) prime component is returned; otherwise all
       summands are returned as a list.

       Warning: If ``prime_decomposition=True`` and
       ``return_all_pieces=False``, then the link returned may have
       fewer components than requested.  This is because a prime piece
       can have fewer components than the link as a whole.


    Some examples:

    >>> K = random_link(25, num_components=1, initial_map_gives_link=True, alternating=True)
    >>> K
    <Link: 1 comp; 25 cross>

    >>> L= random_link(30, consistent_twist_regions=True, simplify = 'global')
    >>> isinstance(random_link(30, return_all_pieces=True), list)
    True
    """
    # This means no trivial loops. PlanarMap accepts 6, which means
    # no bigons, but this is unbearably slow.
    let edge_conn_param = 4

    var link: Link[T]
    if num_components0.isNone:
        link = random_link_internal[T](crossings, edge_conn_param)
    elif initial_map_gives_link:
        link = random_link_internal[T](crossings, edge_conn_param, num_components0.get(), max_tries)
    else:
        let num_components = num_components0.get()
        for i in 0 ..< max_tries:
            link = random_link_internal[T](crossings, edge_conn_param)
            if link.link_components.len >= num_components:
                break
        if link.link_components.len < num_components:
            raise newException(LinkGenerationError, &"did not generate the requested link after {max_tries} tries")
        var comps = longest_components(link, num_components)
        var comps_mask = newSeq[bool](link.link_components.len)
        for comp in comps:
            comps_mask[comp] = true
        link = link.sublink(comps_mask)

    # Adjust the currently random crossings to match the request

    if alternating:
        raise newException(AssertionDefect, "Not implemented")

    if consistent_twist_regions:
        discard make_twist_regions_consistent(link)

    # Initial simplification, if any.

    if simplification_mode.isSome:
        discard link.simplify(simplification_mode.get())

    # Pull into "prime" pieces, if requested.

    if prime_decomposition:
        var prime_pieces = newSeq[Link[T]]()
        for L in link.deconnect_sum():
            var L_copy = L.copy()
            if simplification_mode.isSome:
                discard L_copy.simplify(simplification_mode.get())
            for small_L in L_copy.deconnect_sum():
                prime_pieces.add(small_L)
        if return_all_pieces:
            return prime_pieces
        var largest_prime_piece_index = 0
        for i in 1 ..< prime_pieces.len:
            if prime_pieces[i].crossings.len > prime_pieces[largest_prime_piece_index].crossings.len:
                largest_prime_piece_index = i
        return @[prime_pieces[largest_prime_piece_index]]
    else:
        return @[link]

proc random_link*(
    crossings: int, num_components: Option[int] = none(int),
    initial_map_gives_link: bool = false, alternating: bool = false,
    consistent_twist_regions: bool = false,
    simplification_mode: Option[SimplificationMode] = some(SimplificationMode.basic),
    prime_decomposition: bool = true, return_all_pieces: bool = false,
    max_tries: int = 100
): seq[Link0] =
    return random_link_with_extra_info[int](
        crossings, num_components, initial_map_gives_link, alternating,
        consistent_twist_regions, simplification_mode,
        prime_decomposition, return_all_pieces, max_tries
    )
