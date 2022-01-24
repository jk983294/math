#include <iostream>
#include <discrete_distribution.h>

using namespace std;

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<double> weight = {1.2, 3.4, NAN, 5.6, 7.8, 8.9};
    ornate::discrete_distribution distr(weight);

    for (int j = 0; j < 10; ++j) {
        int choice = distr(gen);
        cout << choice << endl;
    }

    cout << "without placement\n";
    auto ret = ornate::weighted_sample(weight, 4);
    for (size_t j = 0; j < ret.size(); ++j) {
        cout << ret[j] << endl;
    }

    cout << "with placement\n";
    ret = ornate::weighted_sample(weight, 4, true);
    for (size_t j = 0; j < ret.size(); ++j) {
        cout << ret[j] << endl;
    }
}