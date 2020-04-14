////////////////////////////////////////////////////////////////////////
// simulate Radial function in fractional calculus
//
// E(R, pho, zeta, t) = E(R) exp(I(kappa zeta - omega t - m phi))
// E(R)=R^(alfa b) sum((-1)^n lamda^(2 n alfa) d_(2n)R^(alfa n))
// lamda^(2 alfa) = w^2/v_(ph)^2 + (i k )^(2 alfa)
// v_(ph) = 1/sqrt(mu0 epsa0)
// b = +/- m
//
// Abhijit Bhattacharyya
// Dec192019 :: 14:15:57
///////////////////////////////////////////////////////////////////////
#include "fracEM.hh"

// this code just tests Continued Fraction routine
void testContFrac(){
  myComplexD zIn(0.0, 1.0);
  double aIn = 1.8;
  double a0 = 0.0;
  myComplexD cfResult = contFrac(1, a0, aIn, zIn);
  std::cout << " a: " << aIn << "  Z: " << zIn << "   result: " << cfResult << std::endl;
}

// initialize parameters
void initialize(double &freq, double &lambda0, double &omega, int &n4Alpha, double &alpha0, double &alpha1, int &alphaN,
    double &kappa, double &zeta,  double &R0, double &R1, int &rN, int &totN, double &mValue, double &bValue){
    freq = 150.0E+6; // 150 MHz
    lambda0 = vph / freq;
    omega = 2.0 * PI * freq;
    kappa = (2.0 * PI) / lambda0;
    zeta = 1.0; // the z level
    n4Alpha = 2; // 1.0 makes n-2 alpha negative which can't provide gamma function
    alpha0 = 1.0e-3;
    alpha1 = 1.0;
    alphaN = 20;
    R0 = 0.0;
    R1 = 5.0;
    rN = 100;
    totN = 500;
    mValue = 1.425;  // this can be determined later
    bValue = 1.425; // determined from b test function
 }

 void checkArgs(int argc, char **argv){
     if (argc > 3){
         std::cout << " RUN with one argument for testmode 0 or 1" << std::endl;
         std::cout << " RUN with two arguments  testmode=0|1  testM=0|1" << std::endl << std::endl;
         std::cout << " ......................Exiting." << std::endl;
         exit(0);
     }
     if (argc == 1){
         std::cout << " RUN with one argument for testmode 0 or 1" << std::endl;
         std::cout << " RUN with two arguments  testmode=0|1  testM=0|1" << std::endl << std::endl;
         std::cout << " Setting all testmode to 0." << std::endl;
         isTest = 0;
         isTestm = 0;
         return;
     }
     if (argc == 2) {
        std::stringstream ssField(argv[1]);
        int argTest = 0;
        ssField >> argTest;
        isTest = argTest;
        return;
     }
     if (argc == 3) {
         std::stringstream ssField1(argv[1]);
         int argTest1 = 0, argTest2 = 0;
         ssField1 >> argTest1;
         std::stringstream ssField2(argv[2]);
         ssField2 >> argTest2;
         isTest = argTest1;
         isTestm = argTest2;
         return;
     }
 }

 int main(int argc, char** argv){
    checkArgs(argc, argv);
    if (isTestm) mbFileOut.open("mbValues.dat");
    if (!isTestm) rAlphaOut.open("alphaRvalues.dat");

    double freq, omega, mValue;
    double dAlpha, alpha, alpha0, alpha1;
    double R0, R1, delR;
    double kappa, zeta, lambda, lambda0, bVal = 0.0;  // zeta is z level
    myComplexD lambdaP(0.0, 0.0), lambdaTotal(0.0, 0.0);
    int alphaN, rN, n4Alpha, totN;
    double bessJ, bessJD;


    // Function testing blocks ------------------------------------
    //testContFrac(); // testing the conFrac function
    // Function testing block ends ---------------------------------

    // initialize parameters
    initialize(freq, lambda0, omega, n4Alpha, alpha0, alpha1, alphaN, kappa,
        zeta, R0, R1, rN, totN, mValue, bVal);

    // Determine initial parameters after getting initialized variables
    lambda = (x2(omega) / x2(vph)); // myComplexD((x2(omega) / x2(vph)), 0.0);
    dAlpha = std::abs(alpha1 - alpha0) / (double(alphaN));// - 1.0);
    delR = std::abs(R1 - R0) / (double(rN) - 1.0);

    if (isTest) std::cout << " alpha0: " << alpha0 << " alpha1: " << alpha1
    << " alphaN: " << alphaN << " dAlpha: " << dAlpha << std::endl;
    if (isTest) std::cout << " R0: " << R0 << " R1: " << R1 << " rN: " << rN << " delR: "
    << delR << std::endl;
    // start loop over alpha *******
    for (int iAlpha = 0; iAlpha < alphaN; iAlpha++){
      alpha = alpha0 + iAlpha * dAlpha;
      if (isTest) std::cout << " alpha: " << alpha << std::endl;

      if (isTestm) {
          getMbyb(alpha);
          continue;
      }

      lambdaP = lambdaPrime(isTest, n4Alpha, alpha, kappa, zeta);
      lambdaTotal = lambda + lambdaP;
      // start loop over R for each alpha
      for (int iRval = 0; iRval < rN; iRval++){
        double rValue = R0 + iRval * delR;
        double rAlphab = pow(rValue, (alpha * bVal));
        if (isTest) std::cout << " R: " << rValue << " rAlphab: " << rAlphab << std::endl;
        int nPow = std::lround((totN + 1) / 2.0);
        int nVals = 0;
        myComplexD ERSUM(0.0, 0.0);
        for (int nITN = 0; nITN < nPow; nITN++){
            nVals = 2 * nITN;
myBessel(alpha, nVals, bVal, rValue, bessJ, bessJD);

            myComplexD lambda2n = powCMPLX(lambdaTotal, (2.0 * nVals));
            double r2AlphaN = pow(rValue, (2.0 * alpha * nVals));
            double minus1PowN = std::pow(-1.0, nVals);
            myComplexD d2jRecurrence = d2jRecur(isTest, nVals, bVal, alpha, mValue, lambdaTotal);
            ERSUM += (minus1PowN * lambda2n * d2jRecurrence  * r2AlphaN);
            if (isTest) {
                std::cout << "n: " << nVals << " lambda2n: " << lambda2n
                << "  r2AlphaN: " << r2AlphaN << " minus1PowN: " << minus1PowN
                << " d2jRecurrence: " << d2jRecurrence << " ERSUM: " << ERSUM << std::endl;
            }
        } // end of n loop
        ERSUM *= rAlphab;
        std::cout << " alpha: " << alpha << " ERSUM: " << ERSUM << std::endl;
        rAlphaOut << alpha << "   " << rValue << "   " << ERSUM << std::endl;
      }// ****** END of R LOOP
    } // ****** END of alpha LOOP
    return 0;
 }
