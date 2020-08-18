#ifndef __MATH_COMBINATION_H_
#define __MATH_COMBINATION_H_

#include <vector>

namespace ornate {

/**
 * generate next choice,
 * max_choice_count[i] means in position i there are maximum max_choice_count[i] choices available
 * find first pos not equal maximum choice, add 1 and then set backward all 0
 */
inline bool next_choice(std::vector<int>& max_choice_count, std::vector<int>& choice) {
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

}  // namespace ornate

#endif
