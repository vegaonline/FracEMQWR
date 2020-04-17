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
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <complex>
#include <sstream>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/constants/constants.hpp>

typedef std::complex<double> myComplexD;

#define PI        boost::math::constants::pi<double>()
#define mu0       (4.0 * PI * 1.0e-7)
#define epsa0     8.85e-12
#define vph       (1.0 / sqrt(epsa0 * mu0))
#define x2(y)     (y * y)
#define x3(y)     (x2(y) * y)
#define invY(y)   (1.0 / y)
#define xBy2(y)   (0.5 * y)
#define isEven(x) (x % 2 == 0 && x != 0)

bool isTest = 0; // global testing flag
bool isTestm = 0; // get m for series of b

std::ofstream mbFileOut;
std::ofstream rAlphaOut;


// Function declaration
int myFactorial(int);
void checkArgs(int, char**);
myComplexD contFrac(bool, double, double, myComplexD);
void initialize(double &, double &, double &, int &, double &, double &, int &, double &, double &,  double &, double &, int &, int &, double &, double &);
double lambdaPrime(bool, int, double, double, double);
myComplexD powCMPLX(myComplexD, double);
double dFactor(bool, int, double, double, double);
double d2jRecur(bool, int, double, double, double, myComplexD);
void getMbyb(double);
void myBessel(double, int, double, double, double &, double &);

// Function definition


// Determine the Bessel and its derivative
void myBessel(double alpha, int nVal, double bVal, double rVal, double &bessJ, double &bessJD){
    bessJ = 0.0;
    bessJD = 0.0;
    double rBy2 = xBy2(rVal);

    for (int ii = 0; ii < nVal; ii++){
        double numer1 = 0.0, denom1 = 0.0, factor1 = 0.0;
        double numer2 = 0.0, denom2 = 0.0, factor2 = 0.0;
        factor1 = alpha * (2.0 * ii + bVal);
        factor2 = factor1 - 1.0;
        numer1 = std::pow(-1.0, ii) * std::pow(rBy2, factor1);
        denom1 = myFactorial(ii) * (boost::math::tgamma(ii + 1.0 + bVal));
        bessJ += (numer1 / denom1);
        numer2 = std::pow(-1.0, ii) * factor1 * std::pow(rBy2, factor2);
        denom2 = denom1;
        bessJD += (numer2 / denom2);
    }
}

// Determine factorial
int myFactorial(int n) {
    int result = 0;
    result = (n > 1) ? n * myFactorial(n-1) : 1;
    return result;
}

// Determine mPhi using a range of b
void getMbyb(double alpha){
    int numB = 100;
    double bMin = 0.0, bMax = 20.0, delB = 0.0;
    delB = (bMax - bMin) / (numB + 1);
    // int numB = (bMax - bMin) / delB;
    for (int ii = 0; ii < numB; ii++){
        double bValue = bMin + ii * delB;
        double alfab = alpha * bValue;
        double part0 = boost::math::tgamma(alfab + 1.0);
        double part1 = boost::math::tgamma(alpha + 1.0);
        double part2 = x2(part1);
        double part3 = boost::math::tgamma(alfab - 2.0 * alpha + 1.0);
        double part4 = boost::math::tgamma(alfab - alpha + 1.0);
        double m2Value = (part0 / (part2 * part3)) + (part0 / (part4 * part1));
        mbFileOut << 50.0 * alpha << "   " << bValue << "   " << std::sqrt(m2Value) << std::endl;
        // mbFileOut << "   " << bValue << "   " << std::sqrt(m2Value) << alpha << std::endl;
    }
}

// Determine d_2j using recurrence
double d2jRecur(bool isTest, int j, double b, double alpha, double mPhi, double lamTot){
    bool isTest1 = 0;
    double lam2 = x2(lamTot);
    if (j%2) return 0.0; // proceeds for even terms only
    if (j > 0) {
        double d2jRecurRes = d2jRecur(isTest, j - 2, b, alpha, mPhi, lamTot);
        double dFactorRes = dFactor(isTest, j, b, alpha, mPhi);
        // myComplexD result = -1.0 * lam2 * d2jRecurRes * dFactorRes;
        double result = lam2 * d2jRecurRes / dFactorRes; // -1 is taken cared in main
        if (isTest1) std::cout << " d2jRecurResTmp: " << d2jRecurRes << " dFactorRes: " << dFactorRes
            << " d2jRecur_result: " << result
            << " lam2: " << lam2 << " lamTot: " << lamTot
            << std::endl;
        return result;
    } else {
        // return -1.0 * lam2;
        return lam2; // -1 is taken cared in main
    }
}

// Find gamma related denominator in d_2j
double dFactor(bool isTest, int j, double b, double alpha, double mPhi){
    double result = 0.0;
    double part0 = boost::math::tgamma(alpha + 1.0);
    if (isTest) std::cout << " part0: " << part0 << "   ";
    double part1 = boost::math::tgamma(alpha * b +  j * alpha + 1.0);
    if (isTest) std::cout << " part1: " << part1 << "   ";
    double part2 = boost::math::tgamma(alpha * b + (j - 2.0) * alpha + 1.0);
    if (isTest) std::cout << " part2: " << part2 << "   ";
    double part3 = boost::math::tgamma(alpha * b + (j - 1.0) * alpha + 1.0);
    if (isTest) std::cout << " part3: " << part3 << "   ";
    double part4 = x2(mPhi) * x2(part0);
    if (isTest) std::cout << " part4: " << part4 << "   ";
    result = part1 / part2;
    if (isTest) std::cout << " result: " << result << "   ";
    result += ((part0 * part1) / part3);
    if (isTest) std::cout << " result: " << result << "   ";
    result -= part4;
    if (isTest) std::cout << " result: " << result << "   \n";
    return result;
}

// determine Z^x where Z is complex number and x is double
myComplexD powCMPLX(myComplexD c, double num){
    double absVal = std::abs(c);
    double theCoeff1 = pow(x2(absVal), (0.5 * num));
    double theCoeff2 = (std::real(c) == 0.0) ? (0.5 * PI) : atan(std::imag(c) / std::real(c));
    double theRealConv = theCoeff1 * cos(num * theCoeff2);
    double theImagConv = theCoeff1 * sin(num * theCoeff2);
    theRealConv = (std::abs(theRealConv) <= 1.0e-15) ? 0.0 : theRealConv;
    theImagConv = (std::abs(theImagConv) <= 1.0e-15) ? 0.0 : theImagConv;
    myComplexD result(theRealConv, theImagConv);
    return result;
}

// find lambda' part to determine lambda
double lambdaPrime(bool isTest, int n4Alpha, double alpha, double kappa, double zeta){
    double nMinus2Alpha = (double)n4Alpha - 2.0 * alpha;
    myComplexD iKappa = myComplexD(0.0, kappa);
    myComplexD iKappaZeta = iKappa * zeta;
    myComplexD iKappaZetaN = powCMPLX(iKappaZeta, nMinus2Alpha);
    if (isTest) std::cout << " alpha: " << alpha << " ikappaL: " << iKappa
        << " iKappaZetaL: " << iKappaZeta << " (iKappaZeta)^nL: "
        << iKappaZetaN << " n-2*alphaL: " << nMinus2Alpha
        << std::endl;
    myComplexD part1 = powCMPLX(iKappa, (double)n4Alpha) * pow(zeta, nMinus2Alpha);
    if (isTest) std::cout << " part1L: " << part1 << "  ";
    myComplexD part2 = invY(iKappaZetaN);
    if (isTest) std::cout << " part2L: " << part2 << "   ";
    myComplexD part3 = contFrac(isTest, 0.0, nMinus2Alpha, iKappaZeta) * exp(-1.0 * iKappaZeta);
    part3 /= boost::math::tgamma(nMinus2Alpha);
    if (isTest) std::cout << " part3L: " << part3 << "   ";
    double result = std::abs(part1 * (part2 - part3));
    if (isTest) std::cout << " resultL: " << result << std::endl;
    return result;
}

// Continued Fraction
// *** tested working for the a = 0.8, z=(2.0, 0.0) i.e. real number ***
// *** tested also for a=0.8, z=(0.0, 1.0) i.e. comlex number ***
// *** isTest is the testing option 1 to test else 0 ***
myComplexD contFrac(bool isTest, double a0, double alpha, myComplexD zIn){
    int nIter = 6; // total umber of steps for one CF calculation

    // fix ai and bi values for iteration step i
    int itn = nIter;

    myComplexD ai = isEven(itn) ? myComplexD(1.0, 0.0) : zIn;
    double     bi = isEven(itn) ? (0.5 * itn - alpha) : 1.0;
    myComplexD bBYa = bi / ai;

    if (isTest) std::cout << " itn: " << itn << " ai: " << ai << " bi: " << bi << std::endl;

    // start iteration
    for (int itn = (nIter - 1); itn >= 0; itn--){
        ai = isEven(itn) ? myComplexD(1.0, 0.0) : zIn;
        bi = isEven(itn) ? (0.5 * itn - alpha) : 1.0;
        bBYa = bi / (ai + bBYa);
        if (isTest) std::cout << " itn: " << itn << " ai: " << ai
        << " bi: " << bi << " b/a: " << bBYa << std::endl;
    } // End of For loop
    return (a0 + bBYa);
}
