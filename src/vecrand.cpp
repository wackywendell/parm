#include "vecrand.hpp"

engine randengine;
normdistribution mynormaldistribution(0,1);
normgenerator gauss(randengine, mynormaldistribution);
lindistribution mylineardistribution;
lingenerator uniformrand(randengine, mylineardistribution);
// Note: there also exists the 'uniform_on_sphere' distribution
// See http://www.boost.org/doc/libs/1_53_0/doc/html/boost_random/reference.html

// for export, esp. to python, in case you want to be using the same seeds
flt rand01(){return uniformrand();};

#ifdef VEC2D
Vec randVec(){
    return Vec(gauss(), gauss());
}

Vec randVecBoxed(){
    return Vec(uniformrand(), uniformrand());
}
#else
Vec randVec(){
    return Vec(gauss(), gauss(), gauss());
}

Vec randVecBoxed(){
    return Vec(uniformrand(), uniformrand(), uniformrand());
}

Vec randVecSphere(flt radius){
    return randVec().normalized() * (cbrt(uniformrand()) * radius);
}
#endif

unsigned int seed(unsigned int n){
    randengine.seed(n);
    return n;
    
}

unsigned int seed(){
    unsigned int n = static_cast<unsigned int>(time(0));
    randengine.seed(n);
    return n;
}

gaussVec::gaussVec(flt sigma) : distro(0,sigma), gauss(randengine,distro){};


bivariateGauss::bivariateGauss(const flt s1, const flt s2,
        const flt corr) : distro(0,1), gauss(randengine,distro){
            set(s1, s2, corr);
};

void bivariateGauss::set(const flt s1, const flt s2, const flt corr){
    // Taken from Allen and Tildesley, 348
    if(!(s1 >= 0)){
		throw std::invalid_argument("bivariateGauss::set: s1 >= 0");
	} else if(!(s2 >= 0)){
		throw std::invalid_argument("bivariateGauss::set: s2 >= 0");
	} else if(!(corr >= 0)){
		throw std::invalid_argument("bivariateGauss::set: corr >= 0");
	} else if(!(corr <= 1)){
		throw std::invalid_argument("bivariateGauss::set: corr <= 1");
	};
    x11 = s1;
    x21 = s2 * corr;
    x22 = s2 * sqrt(1 - corr*corr);
}

Eigen::Matrix<flt, 1, 2> bivariateGauss::generate(){
    flt x1 = gauss();
    flt x2 = gauss();
    // Taken from Allen and Tildesley, 348
    Eigen::Matrix<flt, 1, 2> p;
    p << x11*x1,
         (x21*x1 + x22*x2);
    return p;
}

VecPair bivariateGauss::genVecs(){
#ifdef VEC2D
    Vec x1 = Vec(gauss(), gauss());
    Vec x2 = Vec(gauss(), gauss());
#else
    Vec x1 = Vec(gauss(), gauss(), gauss());
    Vec x2 = Vec(gauss(), gauss(), gauss());
#endif
    // Taken from Allen and Tildesley, 348
    VecPair p;
    p << x1*x11,
         x1*x21 + x2*x22;
    return p;
}

long double toLD(double e){return (long double) e;};
double fromLD(long double e){return (double) e;};
vector<long double> LDVector(vector<double> dists){
    vector<long double> newdists = vector<long double>();
    for(uint i=0; i<dists.size(); i++){
        newdists.push_back((long double) dists[i]);
    }
    return newdists;
};

#ifdef VEC3D
Matrix best_rotation_matrix(Eigen::Matrix<flt, Eigen::Dynamic, NDIM> &from, Eigen::Matrix<flt, Eigen::Dynamic, NDIM> &to) {
    Eigen::JacobiSVD<Matrix> svd(from.adjoint() * to, Eigen::ComputeFullU | Eigen::ComputeFullV);
    
    Matrix VWprod(svd.matrixV() * svd.matrixU().adjoint());
    if(!VWprod.allFinite()) {
        std::cerr << "BestRotationMatrix ERROR" << std::endl;
        std::cerr << "from:" << std::endl;
        std::cerr << from << std:: endl;
        std::cerr << "to:" << std::endl;
        std::cerr << to << std:: endl;
        
        std::cerr << "U:" << std::endl;
        std::cerr << svd.matrixU() << std:: endl;
        std::cerr << "V:" << std::endl;
        std::cerr << svd.matrixV() << std:: endl;
        std::cerr << "VWprod:" << std::endl;
        std::cerr << VWprod << std:: endl;
    }
    flt det = VWprod.determinant();
    flt d = (det > 0.) ? 1. : 0.;
    
    Vec diagonal_vector;
    for(uint i=0; i<NDIM-1; i++) diagonal_vector(i) = 1.0;
    diagonal_vector(NDIM-1) = d;
    Eigen::DiagonalMatrix<flt, NDIM> diag_d = diagonal_vector.asDiagonal();
    
    Matrix rot = (svd.matrixV()) * diag_d * (svd.matrixU().adjoint());
    return rot;
};
#endif
