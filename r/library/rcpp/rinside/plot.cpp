#include <RInside.h>
#include <unistd.h>

int main(int argc, char *argv[]) {
    // create an embedded R instance
    RInside R(argc, argv);

    // evaluate an R expression with curve()
    // because RInside defaults to interactive=false we use a file
    std::string cmd =
        "tmpf <- tempfile(’curve’); "
        "png(tmpf); "
        "curve(xˆ2, -10, 10, 200); "
        "dev.off();"
        "tmpf";
    // we get the last assignment back, here the filename
    std::string tmpfile = R.parseEval(cmd);

    std::cout << "Could now use plot in " << tmpfile << std::endl;
    unlink(tmpfile.c_str());
    // cleaning up

    // alternatively, by forcing a display we can plot to screen
    cmd = "x11(); curve(xˆ2, -10, 10, 200); Sys.sleep(30);";
    R.parseEvalQ(cmd);
    // parseEvalQ evals without assignment

    return 0;
}
