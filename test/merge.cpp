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
    //arma::cube convCube(probe.n_rows, probe.n_cols, nQ2);
    arma::cube convCube(nQ2, probe.n_rows, probe.n_cols);

    for(int qID = 0; qID < nQ2; ++qID) {
        assert(probe.load(address+name+"_"+to_string(qID) + ".h5", arma::hdf5_binary));
        //convCube.slice(qID) = probe;
        for(int i = 0; i < convCube.n_cols; ++i)
        for(int j = 0; j < convCube.n_slices; ++j)
            convCube(qID, i, j) = probe(i, j);

    }
    convCube.save(address+name+".h5", arma::hdf5_binary);
}



int main()
{
    MergeFiles("convF2");
    MergeFiles("convFL");



    return EXIT_SUCCESS;
}
