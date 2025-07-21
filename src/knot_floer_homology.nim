import std/tables
import std/sugar

import cppstl/std_vector
import cppstl/std_string

{.compile: "ComputeHFKv2/Crossing.cpp".}
{.compile: "ComputeHFKv2/Diagrams.cpp".}
{.compile: "ComputeHFKv2/HomologyRank.cpp".}
{.compile: "ComputeHFKv2/KnotFloer.cpp".}
{.compile: "ComputeHFKv2/Max.cpp".}
{.compile: "ComputeHFKv2/Min.cpp".}
{.compile: "ComputeHFKv2/Report.cpp".}
{.compile: "ComputeHFKv2/Simplify.cpp".}
{.compile: "ComputeHFKv2/Utility.cpp".}
{.compile: "hfk_wrapper.cpp".}

from os import splitPath
{.passC:"-I\"" & currentSourcePath().splitPath.head & "\"" .}

type ReturnedHFK* {.importcpp: "ReturnedHFK", header: "hfk_wrapper.hpp".} = object
    success*: bool
    error_message*: CppString
    modulus*: cint
    ranks_0*: CppVector[cint]
    ranks_1*: CppVector[cint]
    ranks_2*: CppVector[cint]
    total_rank*: cint
    seifert_genus*: cint
    fibered*: bool
    L_space_knot*: bool
    tau*: cint
    nu*: cint
    epsilon*: cint
    generators_0*: CppVector[cint]
    generators_1*: CppVector[cint]
    generators_2*: CppVector[cint]
    differentials_0*: CppVector[cint]
    differentials_1*: CppVector[cint]
    differentials_2*: CppVector[cint]

proc PDCodeToHFKInternal*(pd: CppString, prime: cint): ReturnedHFK {.importcpp: "PDCodeToHFKInternal(@)", header: "hfk_wrapper.hpp".}

type ReturnedMorseCode* {.importcpp: "ReturnedMorseCode", header: "hfk_wrapper.hpp".} = object
    success*: bool
    error_message*: CppString
    events_0*: CppVector[CppString]
    events_1*: CppVector[cint]
    events_2*: CppVector[cint]
    girth*: int

proc PDCodeToMorseInternal*(pd: CppString): ReturnedMorseCode {.importcpp: "PDCodeToMorseInternal(@)", header: "hfk_wrapper.hpp".}

type HFK* = object
    modulus*: int
    ranks*: Table[(int, int), int]
    total_rank*: int
    seifert_genus*: int
    fibered*: bool
    L_space_knot*: bool
    tau*: int
    nu*: int
    epsilon*: int
    generators*: Table[int, (int, int)]
    differentials*: Table[(int, int), int]

proc pd_code_to_hfk*(pd: string, prime: int): HFK =
    let ans0 = PDCodeToHFKInternal(toCppString(pd), cint(prime))
    if not ans0.success:
        raise newException(ValueError, ans0.error_message.toString())
    var ranks = initTable[(int, int), int]()
    var seq_0 = ans0.ranks_0.toSeq()
    var seq_1 = ans0.ranks_1.toSeq()
    var seq_2 = ans0.ranks_2.toSeq()
    for i in 0 ..< seq_0.len:
        ranks[(int(seq_0[i]), int(seq_1[i]))] = int(seq_2[i])
    var generators = initTable[int, (int, int)]()
    seq_0 = ans0.generators_0.toSeq()
    seq_1 = ans0.generators_1.toSeq()
    seq_2 = ans0.generators_2.toSeq()
    for i in 0 ..< seq_0.len:
        generators[int(seq_0[i])] = (int(seq_1[i]), int(seq_2[i]))
    var differentials = initTable[(int, int), int]()
    seq_0 = ans0.differentials_0.toSeq()
    seq_1 = ans0.differentials_1.toSeq()
    seq_2 = ans0.differentials_2.toSeq()
    for i in 0 ..< seq_0.len:
        differentials[(int(seq_0[i]), int(seq_1[i]))] = int(seq_2[i])
    return HFK(modulus: ans0.modulus, ranks: ranks, total_rank: ans0.total_rank,
               seifert_genus: ans0.seifert_genus, fibered: ans0.fibered,
               L_space_knot: ans0.L_space_knot, tau: ans0.tau, nu: ans0.nu, epsilon: ans0.epsilon,
               generators: generators, differentials: differentials)

type MorseCode* = object
    events*: seq[(string, int, int)]
    girth*: int

proc pd_code_to_morse*(pd: string): MorseCode =
    let ans0 = PDCodeToMorseInternal(pd.toCppString())
    if not ans0.success:
        raise newException(ValueError, ans0.error_message.toString())
    var seq_0 = ans0.events_0.toSeq()
    var seq_1 = ans0.events_1.toSeq()
    var seq_2 = ans0.events_2.toSeq()
    var events = collect:
        for i in 0 ..< seq_0.len:
            (seq_0[i].toString(), int(seq_1[i]), int(seq_2[i]))
    return MorseCode(events: events, girth: int(ans0.girth))
