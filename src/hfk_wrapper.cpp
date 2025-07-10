#include <vector>
#include <unordered_map>
#include <string>
#include <mutex>
#include <utility>

#include "ComputeHFKv2/Alg.h"
#include "ComputeHFKv2/Diagrams.h"
#include "hfk_wrapper.hpp"

// Global variables used by functions called from this module.

std::vector<monomial> MonomialStore;
std::unordered_map<monomial, int, Hash> MonomialMap;
int Bridge;
int Modulus;
std::vector<int> UpwardList;
std::vector<int> MatchingList;
std::vector<Arrow> ArrowList;
std::vector<Arrow> NewArrowList;
std::vector<Gen> GeneratorList;
std::vector<Gen> NewGeneratorList;
const monomial MonomialOne = {0};

// Unexported declarations from Alternating.cpp

struct Term
{
    idem Idem;
    int Alexander;
    int Coeff;
};

extern std::vector<Term> AfterMaxAlt(std::vector<Term> Old, int Position);
extern std::vector<Term> AfterMinAlt(std::vector<Term> Old);
extern std::vector<Term> AfterCrossingAlt(std::vector<Term> Old, int Crossing);
extern int Signature (PlanarDiagram Diag);

std::once_flag _monomialStoreAndMapInitialized;

bool isPrime(int n)
{
    if (n < 2)
        return false;
    for (int i = 2; i*i <= n; i++)
        if (n % i == 0)
            return false;
    return true;
}

std::vector<std::pair<std::string, std::pair<int, int>>> MorseListAsEvents(const std::vector<int> &morseList)
{
    std::vector<std::pair<std::string, std::pair<int, int>>> events;

    for (int i = 0; i < sizeAsInt(morseList); i++)
    {
        if (morseList[i] > 999)
        {
            const int c = morseList[++i];
            events.push_back(std::make_pair("cup", make_pair(c - 1, c)));
        }
        else if (morseList[i] > -1000)
        {
            const int c = morseList[i];
            if (c > 0)
                events.push_back(std::make_pair("cross", make_pair(c - 1, c)));
            else
                events.push_back(std::make_pair("cross", make_pair(- c, - c - 1)));
        }
        else
            events.push_back(std::make_pair("cap", make_pair(0, 1)));
    }
    return events;
}

std::pair<std::vector<std::pair<std::string, std::pair<int, int>>>, int> MorseCodeAsEvents(const MorseCode &code)
{
    return std::make_pair(MorseListAsEvents(code.GetMorseList()), code.GetGirth());
}

ReturnedHFK PDCodeToHFKInternal(std::string pd, int prime)
{
    if (prime > 32768)
        return ReturnedHFK("The prime " + std::to_string(prime) + " is greater than 2^15.");
    if (!isPrime(prime))
        return ReturnedHFK(std::to_string(prime) + " is not prime");

    const PlanarDiagram diag = PlanarDiagram(pd);
    if (diag.NotValid())
        return ReturnedHFK("The PD code does not describe a knot projection.");
    if (diag.R1Reducible())
        return ReturnedHFK("The PD code describes a knot projection with a reducing Reidemeister I move");

    const MorseCode lastCheck = diag.GetSmallGirthMorseCode(1);
    if (lastCheck.GetMorseList().size() == 0)
        return ReturnedHFK("This PD code cannot be handled, possibly a connected sum?");

    std::call_once(
        _monomialStoreAndMapInitialized,
        []()
        {
            // Shouldn't we clear this before?
            MonomialStore.push_back(MonomialOne);
            MonomialMap.insert(make_pair(MonomialOne, 0));
        });

    const MorseCode M = diag.GetSmallGirthMorseCode();
    if (M.GetMorseList().empty())
        return ReturnedHFK("Could not compute a small girth Morse code");
    if (M.GetGirth() > 2*MAXBRIDGE)
        return ReturnedHFK("Girth number exceeds " + std::to_string(2 * MAXBRIDGE));

    KnotFloerComplex KFC = ComputingKnotFloer(M, prime, false);
    if (KFC.Prime == 0)
        return ReturnedHFK("Interrupted!");
    return ReturnedHFK(prime,
                       KnotFloerRanks(KFC),
                       KFC.Generators.size(),
                       Genus(KFC),
                       Fibered(KFC),
                       LSpaceKnot(KFC),
                       Tau(KFC),
                       Nu(KFC),
                       Epsilon(KFC),
                       KnotFloerGenerators(KFC),
                       KnotFloerDifferentials(KFC)
                       );
}

ReturnedMorseCode PDCodeToMorseInternal(std::string pd)
{
    const PlanarDiagram diag = PlanarDiagram(pd);
    if (diag.NotValid())
        return ReturnedMorseCode("The PD code does not describe a knot projection.");
    if (diag.R1Reducible())
        return ReturnedMorseCode("The PD code describes a knot projection with a reducing Reidemeister I move");

    const MorseCode lastCheck = diag.GetSmallGirthMorseCode(1);
    if (lastCheck.GetMorseList().size() == 0)
        return ReturnedMorseCode("This PD code cannot be handled, possibly a connected sum?");

    const bool alternating = diag.Alternating();
    // Trying to mimick the behavior of PDCodeToHFK...
    const MorseCode morseCode = alternating ? diag.GetSmallGirthMorseCode(10)
                                            : diag.GetSmallGirthMorseCode();
    if (morseCode.GetMorseList().empty())
        return ReturnedMorseCode("Could not compute a small girth Morse code");
    if (alternating && morseCode.GetGirth() > 2 * MAXBRIDGE)
        return ReturnedMorseCode("Girth number exceeds " + std::to_string(2 * MAXBRIDGE));

    auto return_pair = MorseCodeAsEvents(morseCode);
    auto returned_morse_code = ReturnedMorseCode("");
    returned_morse_code.success = true;
    returned_morse_code.events_0 = std::vector<std::string>();
    returned_morse_code.events_1 = std::vector<int>();
    returned_morse_code.events_2 = std::vector<int>();
    for (const auto& cur_pair : return_pair.first)
    {
        returned_morse_code.events_0.push_back(cur_pair.first);
        returned_morse_code.events_1.push_back(cur_pair.second.first);
        returned_morse_code.events_2.push_back(cur_pair.second.second);
    }
    returned_morse_code.girth = return_pair.second;
    return returned_morse_code;
}
