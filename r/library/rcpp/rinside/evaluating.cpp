#include <RInside.h>

int main(int argc, char *argv[]) {
    RInside R(argc, argv);
    // assignment can be done directly via []
    R["x"] = 10;
    R["y"] = 20;

    // R statement evaluation and result
    R.parseEvalQ("z <- x + y");  // suffix Q means quiet, no return value

    // retrieval access using [] and implicit wrapper
    int sum = R["z"];
    std::cout << "10 + 20 = " << sum << std::endl;

    // we can also return the value directly
    sum = R.parseEval("x + y");
    std::cout << "10 + 20 = " << sum << std::endl;
    return 0;
}
