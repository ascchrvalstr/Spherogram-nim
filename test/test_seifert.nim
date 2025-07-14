import std/unittest

import ../src/links
import ../src/seifert

test "seifert_circles":
    var fig8 = link_from_PD_code(@[[1, 7, 2, 6], [5, 3, 6, 2], [7, 4, 0, 5], [3, 0, 4, 1]])
    check fig8.seifert_circles() == [@[(0, 1), (2, 3), (1, 1), (3, 3)], @[(0, 2), (1, 2)], @[(2, 2), (3, 2)]]
