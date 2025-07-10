#include <vector>
#include <string>
#include <utility>
#include <map>

struct ReturnedHFK
{
    bool success;
    std::string error_message;
    int modulus;
    std::vector<int> ranks_0;
    std::vector<int> ranks_1;
    std::vector<int> ranks_2;
    int total_rank;
    int seifert_genus;
    bool fibered;
    bool L_space_knot;
    int tau;
    int nu;
    int epsilon;
    std::vector<int> generators_0;
    std::vector<int> generators_1;
    std::vector<int> generators_2;
    std::vector<int> differentials_0;
    std::vector<int> differentials_1;
    std::vector<int> differentials_2;

    ReturnedHFK(std::string error_msg) : success{false}, error_message{error_msg} {};

    ReturnedHFK(int modulus0, std::map<std::pair<int, int>, int> ranks0, int total_rank0, int seifert_genus0, bool fibered0, bool L_space_knot0, int tau0, int nu0, int epsilon0, std::map<int, std::pair<int, int>> generators0, std::map<std::pair<int, int>, int> differentials0) : success{true}, modulus{modulus0}, total_rank{total_rank0}, seifert_genus{seifert_genus0}, fibered{fibered0}, L_space_knot{L_space_knot0}, tau{tau0}, nu{nu0}, epsilon{epsilon0}
    {
        ranks_0 = std::vector<int>();
        ranks_1 = std::vector<int>();
        ranks_2 = std::vector<int>();
        for (const auto& cur_pair : ranks0)
        {
            ranks_0.push_back(cur_pair.first.first);
            ranks_1.push_back(cur_pair.first.second);
            ranks_2.push_back(cur_pair.second);
        }
        generators_0 = std::vector<int>();
        generators_1 = std::vector<int>();
        generators_2 = std::vector<int>();
        for (const auto& cur_pair : generators0)
        {
            generators_0.push_back(cur_pair.first);
            generators_1.push_back(cur_pair.second.first);
            generators_2.push_back(cur_pair.second.second);
        }
        differentials_0 = std::vector<int>();
        differentials_1 = std::vector<int>();
        differentials_2 = std::vector<int>();
        for (const auto& cur_pair : differentials0)
        {
            differentials_0.push_back(cur_pair.first.first);
            differentials_1.push_back(cur_pair.first.second);
            differentials_2.push_back(cur_pair.second);
        }
    };
};

ReturnedHFK PDCodeToHFKInternal(std::string pd, int prime);

struct ReturnedMorseCode
{
    bool success;
    std::string error_message;
    std::vector<std::string> events_0;
    std::vector<int> events_1;
    std::vector<int> events_2;
    int girth;

    ReturnedMorseCode(std::string error_message0) : success{false}, error_message{error_message0} {};
};

ReturnedMorseCode PDCodeToMorseInternal(std::string pd);
