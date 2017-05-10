/************************************************************/
//              DESCRIPTION
// Close Fuel cycle scenario :
// This park is constituted by a FBR-Na MOX which
// multi-recycle its own fuel.
// The Storage is initially filled with Pu in order
// to this scenario to be doable
//         _______________	_______     ____    _______
//        |		   |   | FBR   |   |    |  |       |
//  ||===>|FabricationPlant| =>|Reactor| =>|Pool|=>|Storage|===||
//  ||    |________________|   |_______|   |____|  |_______|   ||
//  ||=========================================================||
//
// The spent fuel goes to the pool for 5 y
// then it goes to the Storage
//
//@author BaL
/***********************************************************/
#include <math.h>
#include <iomanip>
#include <sstream>
#include <string>
#include "CLASSHeaders.hxx"

#include "Equivalence/EQM_FBR_BakerRoss_MOX.hxx"
#include "Equivalence/EQM_FBR_MLP_Keff.hxx"
#include "Equivalence/EQM_MLP_Kinf.hxx"
#include "Equivalence/EQM_PWR_MLP_MOX.hxx"
#include "Equivalence/EQM_PWR_MLP_MOX_Am.hxx"
#include "Equivalence/EQM_PWR_POL_UO2.hxx"

using namespace std;

int main(int argc, char** argv) {
  // seconds in one year
  cSecond year = 3600 * 24. * 365.25;
  /******LOG MANAGEMENT**********************************/
  // Definition of the Log file : CLASS messages output
  int Std_output_level = 0;  // Only error are shown in terminal
  int File_output_level =
      2;  // Error + Warning + Info are shown in the file CLASS_OUTPUT.log
  CLASSLogger* Logger =
      new CLASSLogger("CLASS_OUTPUT.log", Std_output_level, File_output_level);

  IsotopicVector IV_1;
  IV_1.Add(92, 230, 0, 2.48825e-16);
  IV_1.Add(92, 231, 0, 3.483e-14);
  IV_1.Add(92, 232, 0, 4.11185e-08);
  IV_1.Add(92, 233, 0, 1.65447e-08);
  IV_1.Add(92, 234, 0, 0.000239692);
  IV_1.Add(92, 235, 0, 2.40811);
  IV_1.Add(92, 236, 0, 0.034307);
  IV_1.Add(92, 237, 0, 0.000210826);
  IV_1.Add(92, 238, 0, 355.988);
  IV_1.Add(93, 237, 0, 0.0944152);
  IV_1.Add(94, 236, 0, 2.23298e-07);
  IV_1.Add(94, 238, 0, 0.0198712);
  IV_1.Add(94, 239, 0, 14.2502);
  IV_1.Add(94, 240, 0, 0.963627);
  IV_1.Add(94, 241, 0, 0.0520374);
  IV_1.Add(94, 242, 0, 0.00225734);
  IV_1.Add(95, 241, 0, 0.018366);
  IV_1.Add(95, 242, 1, 8.8351e-05);
  IV_1.Add(95, 243, 0, 7.68166e-05);
  IV_1.Add(96, 242, 0, 2.63985e-07);
  IV_1.Add(96, 243, 0, 1.6155e-06);
  IV_1.Add(96, 244, 0, 4.92761e-06);
  IV_1.Add(96, 245, 0, 1.98309e-07);
  IV_1.Add(96, 246, 0, 3.7438e-09);
  IV_1 *= 1. / IV_1.GetSumOfAll();

  IsotopicVector IV_2;
  IV_2.Add(92, 235, 0, 3.02497e-05);
  IV_2.Add(92, 238, 0, 0.00417092);
  IV_2 *= 1. / IV_2.GetSumOfAll();

  map<string, IsotopicVector> mymap;
  mymap["Fissile"] = IV_1;
  mymap["Fertile"] = IV_2;

  std::string weight_path =
      "/Users/mouginot/work/MODEL/MODEL_TRU/pwr/EQM/weights/"
      "REP_MOX_TRU_34.2262Wg.xml";
  std::string nfo_file =
      "/Users/mouginot/work/MODEL/MODEL_TRU/pwr/EQM/weights/"
      "REP_MOX_TRU_34.2262Wg.nfo";
  int Batches = 4;
  double k_thresh = 1.034;

  EQM_MLP_Kinf* EQM =
      new EQM_MLP_Kinf(Logger, weight_path, Batches, nfo_file,
                       k_thresh);  // Defining the EquivalenceModel

  double BurnUp = 50;
  // double val = EQM->GetMolarFraction(mymap, BurnUp)["Fissile"];

  IsotopicVector IV_test;
  IV_test.Add(92, 234, 0, 0.00000E+00);
  IV_test.Add(92, 235, 0, 4.29601E+00);
  IV_test.Add(92, 236, 0, 0.00000E+00);
  IV_test.Add(92, 238, 0, 5.99925E+02);
  IV_test.Add(93, 237, 0, 2.14079E-01);
  IV_test.Add(94, 236, 0, 5.04170E-07);
  IV_test.Add(94, 238, 0, 4.52466E-02);
  IV_test.Add(94, 239, 0, 3.25843E+01);
  IV_test.Add(94, 240, 0, 2.21265E+00);
  IV_test.Add(94, 241, 0, 1.19986E-01);
  IV_test.Add(94, 242, 0, 5.22653E-03);
  IV_test.Add(95, 241, 0, 4.23477E-02);
  IV_test.Add(95, 242, 1, 2.04564E-04);
  IV_test.Add(95, 243, 0, 1.78594E-04);
  IV_test.Add(96, 242, 0, 6.11216E-07);
  IV_test.Add(96, 243, 0, 3.75593E-06);
  IV_test.Add(96, 244, 0, 1.15036E-05);
  IV_test.Add(96, 245, 0, 4.64857E-07);
  IV_test.Add(96, 246, 0, 8.81176E-09);
  IV_test *= 1 / IV_test.GetSumOfAll();

  IsotopicVector IV_test_2;
  IV_test_2.Add(92, 235, 0, 2.53029);
  IV_test_2.Add(92, 238, 0, 348.884);
  IV_test_2.Add(93, 237, 0, 0.137248);
  IV_test_2.Add(94, 236, 0, 3.246e-07);
  IV_test_2.Add(94, 238, 0, 0.028886);
  IV_test_2.Add(94, 239, 0, 20.7149);
  IV_test_2.Add(94, 240, 0, 1.40079);
  IV_test_2.Add(94, 241, 0, 0.0756447);
  IV_test_2.Add(94, 242, 0, 0.00328141);
  IV_test_2.Add(95, 241, 0, 0.026698);
  IV_test_2.Add(95, 242, 1, 0.000128432);
  IV_test_2.Add(95, 243, 0, 0.000111665);
  IV_test_2.Add(96, 242, 0, 3.83745e-07);
  IV_test_2.Add(96, 243, 0, 2.34839e-06);
  IV_test_2.Add(96, 244, 0, 7.16308e-06);
  IV_test_2.Add(96, 245, 0, 2.88274e-07);
  IV_test_2.Add(96, 246, 0, 5.44222e-09);
  IV_test_2 *= 1 / IV_test_2.GetSumOfAll();

  double BU_t = EQM->GetMaximumBurnUp_MLP(IV_test, 50);
  std::cout << "bu-test " << BU_t << std::endl;
  double BU_2 = EQM->GetMaximumBurnUp_MLP(IV_test_2, 50);
  std::cout << "bu-2 " << BU_2 << std::endl;

  double BU_1 = EQM->GetMaximumBurnUp_MLP(IV_1, 50);
  std::cout << "bu-1 " << BU_1 << std::endl;
  IV_test.Print();
  IV_1.Print();
  (IV_test - IV_1).Print();
}

//==========================================================================================
// Compilation
//==========================================================================================
/*

 \rm CLASS* ; g++ -o CLASS_Exec model_test.cxx -I $CLASS_include -L $CLASS_lib
 -lCLASSpkg `root-config --cflags` `root-config --libs` -fopenmp -lgomp
 -Wunused-result


 */
