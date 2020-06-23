#include <RInside.h>  // for the embedded R via RInside
#include <mpi.h>      // mpi header

int main(int argc, char *argv[]) {
    // mpi initialization
    MPI::Init(argc, argv);

    // obtain current node rank and total nodes running
    int myrank = MPI::COMM_WORLD.Get_rank();
    int nodesize = MPI::COMM_WORLD.Get_size();

    // create an embedded R instance
    RInside R(argc, argv);

    std::stringstream txt;
    // node information
    txt << "Hello from node " << myrank << " of " << nodesize << " nodes!" << std::endl;
    // assign string var to R variable ’txt’
    R.assign(txt.str(), "txt");

    // show node information
    std::string evalstr = "cat(txt)";
    // eval the init string, ignoring any returns
    R.parseEvalQ(evalstr);

    // mpi finalization
    MPI::Finalize();
    return 0;
}
