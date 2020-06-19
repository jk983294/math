#include <L1TrendFilter.h>

using namespace ornate;

int main(int argc, char *argv[]) {
    const char *ifile_y = "/home/aaron/github/math/data/snp500.txt";

    FILE *fp = fopen(ifile_y, "r");
    if (fp == nullptr) {
        printf("ERROR: Could not open file: %s\n", ifile_y);
        exit(EXIT_FAILURE);
    }

    std::vector<double> data;
    double val;
    while (fscanf(fp, "%lg\n", &val) != EOF) {
        data.push_back(val);
    }
    fclose(fp);

    ornate::L1TrendFilter filter;
    filter.set_data(data);
    filter.calc();

    printf("index,y,x\n");
    for (size_t i = 0; i < filter.m_y.size(); ++i) {
        printf("%zu,%.6f,%.6f\n", i, filter.m_y[i], filter.m_x[i]);
    }
}
