#include <armadillo>
#include <string>
#include <cassert>

using namespace std;

void MergeFiles(string name)
{
    string address = "/nfs/dust/cms/user/zlebcr/Krakow/convMat/";
    //string name= "hoho";
    int nQ2 = 46;

    cout << "Loading file " << address+name+"_0.h5" << endl;
    arma::mat probe;
    assert(probe.load(address+name+"_0.h5", arma::hdf5_binary));

    cout << "N rows is " << probe.n_rows << endl;
    cout << "N cols is " << probe.n_cols << endl;
    arma::cube convCube(probe.n_rows, probe.n_cols, nQ2);

    for(int i = 0; i < nQ2; ++i) {
        assert(probe.load(address+name+"_"+to_string(i) + ".h5", arma::hdf5_binary));
        convCube.slice(i) = probe;

    }
    convCube.save(address+name+".h5", arma::hdf5_binary);
}



int main()
{
    MergeFiles("convFT");
    MergeFiles("convFL");



    return EXIT_SUCCESS;
}
