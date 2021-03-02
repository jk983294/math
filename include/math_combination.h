#ifndef __MATH_COMBINATION_H_
#define __MATH_COMBINATION_H_

#include <vector>

namespace ornate {

/**
 * generate next choice,
 * max_choice_count[i] means in position i there are maximum max_choice_count[i] choices available
 * find first pos not equal maximum choice, add 1 and then set backward all 0
 */
inline bool next_choice(const std::vector<int>& max_choice_count, std::vector<int>& choice) {
    int size = static_cast<int>(max_choice_count.size());
    for (int i = size - 1; i >= 0; --i) {
        if (choice[i] != max_choice_count[i] - 1) {
            ++choice[i];
            fill(choice.begin() + i + 1, choice.end(), 0);
            return true;
        }
    }
    return false;
}

inline std::vector<std::vector<int>> get_all_choices(const std::vector<int>& max_choice_count) {
    std::vector<std::vector<int>> all_choices;
    std::vector<int> choice(max_choice_count.size(), 0);
    do {
        all_choices.push_back(choice);
    } while (ornate::next_choice(max_choice_count, choice));
    return all_choices;
}

inline std::vector<std::vector<int>> get_all_choices(int n, int max_choice_cnt) {
    std::vector<int> max_choice_count(n, max_choice_cnt);
    return get_all_choices(max_choice_count);
}

}  // namespace ornate

#endif
