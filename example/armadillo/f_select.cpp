#include <math_feature.h>
#include <math_random.h>
#include <vector>

using namespace std;

int main() {
    size_t n = 10;
    ornate::ForwardStep fs;
    fs.set_remove_na_ret(true);
    fs.set_has_intercept(true);
    // fs.set_criterion("aic");
    vector<double> ret = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f0 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f1 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f2 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f3 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f4 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f5 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
//    printf("ret:%f,%f,%f,%f,%f\n", ret[0], ret[1], ret[2], ret[3], ret[4]);
//    printf("%f,%f,%f,%f,%f\n", f0[0], f0[1], f0[2], f0[3], f0[4]);
//    printf("%f,%f,%f,%f,%f\n", f1[0], f1[1], f1[2], f1[3], f1[4]);
//    printf("%f,%f,%f,%f,%f\n", f2[0], f2[1], f2[2], f2[3], f2[4]);
    fs.set_y(ret.data(), n);
    fs.add_data(0, f0.data());
    fs.add_data(1, f1.data());
    fs.add_data(2, f2.data());
    fs.add_data(3, f3.data());
    fs.add_data(4, f4.data());
    fs.add_data(5, f5.data());
    fs.process();
    auto selected = fs.get_selected();
    printf("selected=");
    for (int s : selected) {
        printf("%d,", s);
    }
    printf("\n");
    return 0;
}
