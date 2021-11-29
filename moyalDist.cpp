#include "moyalDist.hh"

double moyalDist::sampleInvCDF(double mean, double sigma, double rand) const{

    // The Moyal distribution cdf is cdf(x; mu, sigma) = erfc(exp(- 1/2 * (x - mu) / sigma) / sqrt(2))
    // This is an implementation of erfc^-1 (z) from https://functions.wolfram.com/06.31.06.0002.01
    // kmax is the number of terms in the series to keep. Higher the kmax better is the agreement with accept/reject m
    // but high kmax also incurs high computational cost

    const int kmax = 20; 
    //std::vector<double> coeff;
    coeff.push_back(1); // coeff[0] = 1; 
    double t = 0.5 * std::sqrt(std::acos(-1.)) * (1 - rand);
    double sum = 0;

    for(int k=0; k < kmax; k++ ){

        double sum_coeff = 0;
        if( k > 0 ){
            for(int m = 0; m < k; m++){
                sum_coeff += ( coeff.at(m) * coeff.at( k - 1 - m) ) / ( (m + 1) * (2.0 * m + 1) );
            }
            coeff.push_back(sum_coeff); // these coeffcients can be pre-calculated and used once kmax is fixed
             
        }
        sum +=  ( coeff.at(k) / (2.* k + 1.0) ) * std::pow(t, (2.* k + 1.0));
        // std::cout << "k = " << k << " && c = " <<  coeff.at(k) / (2.* k + 1.0) << std::endl;
    } 
    
    double y = sum;
    double x = mu - 2 * sigma * std::log ((std::sqrt(2.) * y));
    return x;
}


double moyalDist::sampleAR(double mean, double sigma) const{

    // This is an implementation of accep-reject method for sampling Moyal distribution
    // recipe at http://www.stat.rice.edu/~dobelman/textfiles/DistributionsHandbook.pdf
    // Needs two random numbers. 

    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);
    const double pi = std::acos(-1);
    const double hmax = 0.912;
    double h = 0; 
    double Hy = 0;
    double z = 0;

    do
    {
        double rand1 = dis(gen);
        double rand2 = dis(gen);
        //std::cout << "rand1 == " << rand1 << " && " << "rand2 == " << rand2 << "\n";

        double y = pi * (rand1 - 0.5);
        h = rand2 * hmax;
        z = std::tan(y);
        //std::cout << "z == " << z << "\n" ;
        Hy = std::sqrt(1.0/(2.0 * pi)) * (1.0 + std::pow(z, 2) ) \
                    * std::exp( -0.5 * ( z + std::exp(-z) ) );

    } while ( h > Hy );
    double x = mu + sigma * z;
    return (x) ;   
}