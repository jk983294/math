#include <math_feature.h>
#include <thread>
#include <vector>

using namespace std;

void generate_random_data(vector<double>& data, size_t n, mt19937& generator) {
    uniform_real_distribution<double> uid(-1.0, 1.0);
    data.resize(n);
    for (size_t i = 0; i < n; ++i) {
        data[i] = uid(generator);
    }
}

int main() {
    random_device rd;  // non-deterministic generator
    mt19937 generator(42);
    size_t n = 10000000;
    int f_cnt = 10;
    vector<vector<double>> features;

    ornate::ForwardStep fs;
    fs.set_thread_num(std::thread::hardware_concurrency());
    fs.set_remove_na_ret(true);
    fs.set_has_intercept(true);
    // fs.set_criterion("aic");
    features.push_back(vector<double>(n, 0));
    auto& ret = features.back();
    generate_random_data(ret, n, generator);
    fs.set_y(ret.data(), n);

    for (int i = 0; i < f_cnt; ++i) {
        features.push_back(vector<double>(n, 0));
        auto& feature = features.back();
        generate_random_data(feature, n, generator);
        fs.add_data(i, feature.data());
    }

    fs.process();
    auto selected = fs.get_selected();
    printf("selected=");
    for (int s : selected) {
        printf("%d,", s);
    }
    printf("\n");
    return 0;
}
