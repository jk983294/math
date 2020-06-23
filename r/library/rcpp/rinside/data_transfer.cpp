#include <RInside.h>

int main(int argc, char *argv[]) {
    RInside R(argc, argv);
    double d1 = 1.234;

    std::vector<double> d2;
    d2.push_back(1.23);
    d2.push_back(4.56);

    std::map<std::string, double> d3;  // map
    d3["a"] = 7.89;
    d3["b"] = 7.07;

    std::list<double> d4;  // list of doubles
    d4.push_back(1.11);
    d4.push_back(4.44);

    R["d1"] = d1;  // scalar double
    R["d2"] = d2;  // vector of doubles
    R["d3"] = d3;
    R["d4"] = d4;

    // now access in R
    std::string txt =
        "cat('\nd1=', d1, '\n'); print(class(d1));"
        "cat('\nd2=\n'); print(d2); print(class(d2));"
        "cat('\nd3=\n'); print(d3); print(class(d3));"
        "cat('\nd4=\n'); print(d4); print(class(d4));";

    R.parseEvalQ(txt);
    return 0;
}
