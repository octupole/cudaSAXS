#include "readTraj.h"
bool readTraj::firsttime = true;
readTraj::readTraj(string filename, int numAtoms)
{
    // Open the XTC file
    matrix box;

    fin = xdrfile_open(filename.c_str(), "r");
    fp = xdrfile_get_fp(fin);
    xdr = xdrfile_get_xdr(fin);
    int bOK = 0;
    // Read the number of atoms
    if (read_xtc_natoms(filename.c_str(), &natoms) != exdrOK)
    {
        std::cerr << "Error: Unable to read number of atoms from " << filename << std::endl;
        return;
    }
    x0 = new rvec[natoms];
    coords = vector<vector<float>>(natoms, vector<float>(3));
    int next = xtc_get_next_frame_number(fp, xdr, natoms);
    dstep = next;

    cout << dstep << endl;
}
void readTraj::readaStep(int frameO)
{
    int step;
    float timeA;
    float prec;
    matrix box;

    int yy{0};
    if (firsttime)
    {
        firsttime = false;
    }
    int frame = frameO * dstep;

    do
    {
        int status = read_xtc(fin, natoms, &step, &timeA, box, x0, &prec);
        if (status == exdrENDOFFILE)
        {
            std::cerr << "Error: Reached end of file before finding desired frame." << std::endl;
            break;
        }
        else if (status != exdrOK)
        {
            std::cerr << "Error: Failed to read frame." << std::endl;
            break;
        }
    } while (frame != step);
    cout << "Frame: " << frame << " Step: " << step << " Time: " << timeA << " Prec: " << prec << endl;
    step_c = step;
    time_c = timeA;

    setCoord(x0);
}
// try
// {
//     if (firsttime)
//     {
//         firsttime = false;
//         if (my_nframe != fin->goffStep())
//         {
//             if (!(my_nframe <= fin->gFrameNumber()))
//                 throw string(" Beginning frame is out of range ");
//         }
//         fin->Rewind();
//     }
//     // Apparenly xdr_xtc_seek_frame does not work for the beginning frame!!
//     if (my_nframe != fin->goffStep())
//         xdr_xtc_seek_frame(my_nframe, fp, myxdr, natom);
//     if (read_xtc(xd, natom, &step, &time, box, x0, &prec))
//         throw string("Cannot read next frame!");
// }
// catch (const string &s)
// {
//     cout << s << endl;
//     exit(-1);
// }
void readTraj::setCoord(rvec *x0)
{

    for (auto i = 0; i < natoms; i++)
        for (auto j = 0; j < DIM; j++)
            coords[i][j] = x0[i][j] * Units::nanoToA;
}
vector<vector<float>> readTraj::getInverseMatrix(vector<vector<float>> &co)
{
    vector<vector<float>> T_inv{3, vector<float>(3)};

    // Check for zero diagonal elements
    if (co[0][0] == 0 || co[1][1] == 0 || co[2][2] == 0)
    {
        throw std::runtime_error("Matrix is singular and cannot be inverted");
    }

    // Calculate the inverse elements
    T_inv[0][0] = 1.0 / co[0][0];
    T_inv[0][1] = -co[0][1] / (co[0][0] * co[1][1]);
    T_inv[0][2] = (co[0][1] * co[1][2] - co[0][2] * co[1][1]) / (co[0][0] * co[1][1] * co[2][2]);

    T_inv[1][1] = 1.0 / co[1][1];
    T_inv[1][2] = -co[1][2] / (co[1][1] * co[2][2]);

    T_inv[2][2] = 1.0 / co[2][2];

    return T_inv;
}
readTraj::~readTraj()
{
    delete[] x0;
}