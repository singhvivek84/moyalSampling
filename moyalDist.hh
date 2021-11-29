#include <vector>

class moyalDist{
    public:
        double getMean(double mpv, double xi ) const;
        double getSigma(double mpv, double xi) const;
        double sampleMoyalAR(double mpv, double xi) const;
        double sampleMoyalCDF(double mpv, double xi, double rand) const;
        double sampleMoyalCDF(double mpv, double xi, double rand, const int max) const;
    private:
        double mean;
        double sigma; 
        double kmax = 20;
        std::vector<double> coeff;   
};