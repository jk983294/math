#include <math_random.h>
#include <iostream>

using namespace ornate;
using namespace std;

int main() {
    vector<int> ret = ornate::sample_by_n(10, 3);
    printf("%d,%d,%d\n", ret[0], ret[1], ret[2]);
    ret = ornate::sample_by_n(10, 3);
    printf("%d,%d,%d\n", ret[0], ret[1], ret[2]);

    ret = ornate::sample_by_ratio(10, 0.3);
    printf("%d,%d,%d\n", ret[0], ret[1], ret[2]);
    ret = ornate::sample_by_ratio(10, 0.3);
    printf("%d,%d,%d\n", ret[0], ret[1], ret[2]);

    ret = ornate::sample_by_n(10, 20);
    printf("%d,%d,%zu\n", ret[0], ret[9], ret.size());
    return 0;
}
