#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <fstream>
#include <eigen3/Eigen/Dense>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <boost/multi_array.hpp>

using namespace std;
using namespace boost;
using namespace Eigen;
using boost::format;
using boost::extents;
namespace po = boost::program_options;

typedef Array<short int, Dynamic, Dynamic> ArrayXXs;

typedef boost::multi_array<short int, 3> ArrayXXXs;

namespace Util {
    template <typename T>
    inline T get(Array<T, Dynamic, Dynamic> arr, int row, int col) {
        int rows = arr.rows(), cols = arr.cols();
        return arr((row+rows)%rows, (col+cols)%cols);
    }

    template <typename T>
    inline T get(boost::multi_array<T, 3> arr, int row, int col, int dep, int N) {
        return arr[(row+N)%N][(col+N)%N][(dep+N)%N];
    }

    void plusMinus(int N, ArrayXXi &mat, const ArrayXXXs &sigma) {
        if (sigma[0][0][0] > 0)
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    mat(i, j) += sigma[0][i][j];
        else
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j)
                    mat(i, j) -= sigma[0][i][j];
    }
}

class RandPair {
private:
    std::mt19937 gen;
    std::uniform_int_distribution<> dis;
public:
    RandPair(const int l, const int r) {
        std::random_device rd;
        gen = std::mt19937(rd());
        dis = std::uniform_int_distribution<>(l, r);
    }

    int generate() {
        return dis(gen);
    }
};

struct StateTuple {
    int intersum, sum;
    StateTuple(const int intersum, const int sum): intersum(intersum), sum(sum) {}
};

class StateBase {
public:
    virtual const double getInterSum() const = 0;
    virtual const double getSum() const = 0;
    virtual const int getProdNum() const = 0;
    virtual const StateTuple gibbsSampling(int rowId, int colId, int depId, int groupId) const = 0;
    virtual void flip(int rowId, int colId, int depId, int groupId, const StateTuple &tuple) = 0;
    virtual void updateCorr(ArrayXXi &mat_ii, ArrayXXi &mat_hh) const = 0;
    virtual void print(std::ostream &out) const = 0;
};

class State : public StateBase{
private:
    ArrayXXs sigma_ii, sigma_hh;
    int intersum, sum;

    void setInterSum() {
        int cols = sigma_ii.cols(), rows = sigma_ii.rows();

        intersum = 0;

        intersum += (sigma_ii * sigma_hh).sum();

        intersum += (sigma_ii.rightCols(cols-1) * sigma_hh.leftCols(cols-1)).sum();
        intersum += (sigma_ii.col(0) * sigma_hh.col(cols-1)).sum();

        intersum += (sigma_ii.bottomRows(rows-1) * sigma_hh.topRows(rows-1)).sum();
        intersum += (sigma_ii.row(0) * sigma_hh.row(rows-1)).sum();

        intersum += (sigma_ii.bottomRightCorner(rows-1, cols-1) * sigma_hh.topLeftCorner(rows-1, cols-1)).sum();
        intersum += (sigma_ii.bottomLeftCorner(rows-1, 1) * sigma_hh.topRightCorner(rows-1, 1)).sum();
        intersum += (sigma_ii.topRightCorner(1, cols-1) * sigma_hh.bottomLeftCorner(1, cols-1)).sum();
        intersum += sigma_ii(0, 0) * sigma_hh(rows-1, cols-1);
    }

    void setSum() {
        sum = sigma_ii.sum() + sigma_hh.sum();
    }
    
    State() {}
public:
    State(const int N) {
        sigma_ii = ArrayXXs::Constant(N, N, 0);
        sigma_ii.setRandom();
        sigma_ii = (2*sigma_ii+1).sign();
        sigma_hh = ArrayXXs::Constant(N, N, 0);
        sigma_hh.setRandom();
        sigma_hh = (2*sigma_hh+1).sign();
        
        this->setInterSum();
        this->setSum();
    }

    const double getInterSum() const { return intersum; }
    
    const double getSum() const { return sum; }
    
    const int getProdNum() const { return sigma_ii.size() * 2; }

    const StateTuple gibbsSampling(int rowId, int colId, int depth, int groupId) const {
        int type = groupId % 2;
        int s, intersum_;
        switch (type) {
            case 0: // ii
                s = -sigma_ii(rowId, colId);
                intersum_ = intersum + 2*s*(
                    Util::get(sigma_hh, rowId  , colId-1)
                  + Util::get(sigma_hh, rowId-1, colId-1)
                  + Util::get(sigma_hh, rowId  , colId  )
                  + Util::get(sigma_hh, rowId-1, colId  )
                );
                break;
            case 1: // hh
                s = -sigma_hh(rowId, colId);
                intersum_ = intersum + 2*s*(
                    Util::get(sigma_ii, rowId  , colId+1)
                  + Util::get(sigma_ii, rowId+1, colId+1)
                  + Util::get(sigma_ii, rowId  , colId  )
                  + Util::get(sigma_ii, rowId+1, colId  )
                );
                break;
        }
        int sum_ = sum + 2*s;

        return StateTuple(intersum_, sum_);
    }

    void flip(int rowId, int colId, int depth, int groupId, const StateTuple &tuple) {
        int type = groupId % 2;
        switch (type) {
            case 0: // ii
                sigma_ii(rowId, colId) *= -1;
                break;
            case 1: // hh
                sigma_hh(rowId, colId) *= -1;
                break;
        }
        intersum = tuple.intersum;
        sum = tuple.sum;
    }

    void updateCorr(ArrayXXi &mat_ii, ArrayXXi &mat_hh) const {
        if (sigma_ii(0, 0) > 0) {
            mat_ii += sigma_ii.cast<int>();
            mat_hh += sigma_hh.cast<int>();
        } else {
            mat_ii -= sigma_ii.cast<int>();
            mat_hh -= sigma_hh.cast<int>();
        }
    }

    void print(std::ostream &out) const {
        out
            << sigma_ii
            << endl
            << sigma_hh
            << endl
            ;
    }

    friend std::ostream& operator<<(std::ostream &out, const State &s) {
        out
            << format("State @(intersum=%f, sum=%f)") % s.intersum % s.sum
            << endl
            ;
        s.print(out);
        return out;
    }
};

class State3d : public StateBase{
private:
    ArrayXXXs sigma_ii, sigma_hh;
    int N;
    int intersum, sum;
    
    State3d() {}
public:
    State3d(const int N): N(N) {
        sigma_ii.resize(extents[N][N][N]);
        sigma_hh.resize(extents[N][N][N]);
        std::fill_n(sigma_ii.data(), sigma_ii.num_elements(), 1);
        std::fill_n(sigma_hh.data(), sigma_hh.num_elements(), 1);
        this->intersum = 8*N*N*N;
        this->sum = 2*N*N*N;
    }

    const double getInterSum() const { return intersum; }
    
    const double getSum() const { return sum; }
    
    const int getProdNum() const { return sigma_ii.num_elements() * 2; }

    const StateTuple gibbsSampling(int rowId, int colId, int depId, int groupId) const {
        int type = groupId % 2;
        int s, intersum_;
        switch (type) {
            case 0: // ii
                s = -sigma_ii[rowId][colId][depId];
                intersum_ = intersum + 2*s*(
                    Util::get(sigma_hh, rowId-1, colId-1, depId-1, N)
                  + Util::get(sigma_hh, rowId  , colId-1, depId-1, N)
                  + Util::get(sigma_hh, rowId-1, colId  , depId-1, N)
                  + Util::get(sigma_hh, rowId  , colId  , depId-1, N)
                  + Util::get(sigma_hh, rowId-1, colId-1, depId  , N)
                  + Util::get(sigma_hh, rowId  , colId-1, depId  , N)
                  + Util::get(sigma_hh, rowId-1, colId  , depId  , N)
                  + Util::get(sigma_hh, rowId  , colId  , depId  , N)
                );
                break;
            case 1: // hh
                s = -sigma_hh[rowId][colId][depId];
                intersum_ = intersum + 2*s*(
                    Util::get(sigma_ii, rowId+1, colId+1, depId+1, N)
                  + Util::get(sigma_ii, rowId  , colId+1, depId+1, N)
                  + Util::get(sigma_ii, rowId+1, colId  , depId+1, N)
                  + Util::get(sigma_ii, rowId  , colId  , depId+1, N)
                  + Util::get(sigma_ii, rowId+1, colId+1, depId  , N)
                  + Util::get(sigma_ii, rowId  , colId+1, depId  , N)
                  + Util::get(sigma_ii, rowId+1, colId  , depId  , N)
                  + Util::get(sigma_ii, rowId  , colId  , depId  , N)
                );
                break;
        }
        int sum_ = sum + 2*s;

        return StateTuple(intersum_, sum_);
    }

    void flip(int rowId, int colId, int depId, int groupId, const StateTuple &tuple) {
        int type = groupId % 2;
        switch (type) {
            case 0: // ii
                sigma_ii[rowId][colId][depId] *= -1;
                break;
            case 1: // hh
                sigma_hh[rowId][colId][depId] *= -1;
                break;
        }
        intersum = tuple.intersum;
        sum = tuple.sum;
    }

    void updateCorr(ArrayXXi &mat_ii, ArrayXXi &mat_hh) const {
        // Util::plusMinus(N, mat_ii, sigma);
    }

    void print(std::ostream &out) const {
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                out << sigma_ii[i][j][0] << " ";
            }
            out << endl;
        }
        out << endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                out << sigma_hh[i][j][0] << " ";
            }
            out << endl;
        }
        out << endl;
    }

    friend std::ostream& operator<<(std::ostream &out, const State3d &s) {
        out
            << format("State3d @(intersum=%f, sum=%f)") % s.intersum % s.sum
            << endl
            ;
        s.print(out);
        return out;
    }
};

class System {
private:
    double J, h;
public:
    System(const double J, const double h): J(J), h(h) {}

    double hamiltonian(const StateBase &state) const {
        double intersum = state.getInterSum(), sum = state.getSum();
        return -(J*intersum + h*sum);
    }

    double hamiltonian(const StateTuple &tuple) const {
        double intersum = tuple.intersum, sum = tuple.sum;
        return -(J*intersum + h*sum);
    }
};

class Measurement {
private:
    double energy, magnetization;
public:
    Measurement(const System &system, const StateBase &state, const double beta) {
        energy = system.hamiltonian(state);
        magnetization = state.getSum() / state.getProdNum();
    }
};

struct Option {
    int outputInt = 1500, skip = 200, maxItr = 1e5;
    double J = 1.0, h = 0.0, beta = 1.0;
    int N = 20, dimension = 2;
    string corrDump = "", confDump = "";
};

Option argument_handler(int argc, char *argv[]) {
    // Declare the supported options.
    po::options_description desc("Options are");
    desc.add_options()
        ("outputInt", po::value<int>(), "output interval between successive uncorrelated samples")
        ("skip", po::value<int>(), "skipped uncorrelated samples")
        ("maxItr", po::value<int>(), "max iteration limit")
        ("J", po::value<double>(), "J")
        ("h", po::value<double>(), "h")
        ("beta", po::value<double>(), "beta")
        ("N", po::value<int>(), "system size")
        ("dimension", po::value<int>(), "system dimension")
        ("corrDump", po::value<string>(), "correlation dump file")
        ("confDump", po::value<string>(), "configuration dump file")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        exit(1);
    }

    Option option;
    if (vm.count("outputInt")) option.outputInt = vm["outputInt"].as<int>   ();
    if (vm.count("skip"))      option.skip      = vm["skip"]     .as<int>   ();
    if (vm.count("maxItr"))    option.maxItr    = vm["maxItr"]   .as<int>   ();
    if (vm.count("J"))         option.J         = vm["J"]        .as<double>();
    if (vm.count("h"))         option.h         = vm["h"]        .as<double>();
    if (vm.count("beta"))      option.beta      = vm["beta"]     .as<double>();
    if (vm.count("N"))         option.N         = vm["N"]        .as<int>   ();
    if (vm.count("dimension")) option.dimension = vm["dimension"].as<int>   ();
    if (vm.count("corrDump"))  option.corrDump  = vm["corrDump"] .as<string>();
    if (vm.count("confDump"))  option.confDump  = vm["confDump"] .as<string>();

    cout << "Parameter list:\n"
         << "outputInt: " << option.outputInt << endl
         << "skip:      " << option.skip      << endl
         << "maxItr:    " << option.maxItr    << endl
         << "J:         " << option.J         << endl
         << "h:         " << std::setprecision(5) << option.h         << endl
         << "beta:      " << option.beta      << endl
         << "N:         " << option.N         << endl
         << "dimension: " << option.dimension << endl
         << "confDump:  " << option.confDump  << endl
         << "corrDump:  " << option.corrDump  << endl;

    return option;
}

int main(int argc, char *argv[]) {
    const Option option = argument_handler(argc, argv);

    // Init
    StateBase *state;
    if (option.dimension == 2) {
        state = new State(option.N);
    } else {
        state = new State3d(option.N);
    }
    const System system(option.J, option.h);
    double hamiltonian = system.hamiltonian(*state);
    ArrayXXi corrMat_ii = ArrayXXi::Constant(option.N, option.N, 0),
        corrMat_hh = ArrayXXi::Constant(option.N, option.N, 0);

    int itr = 0;
    const int coutInt = option.maxItr/100;
    RandPair randPair(0, option.N-1);
    RandPair randPair_g(0, 24-1);
    Array4i randPlacer;

    std::random_device unif_rd;
    std::mt19937 unif_gen(unif_rd());
    std::uniform_real_distribution<> unif_dis(0.0, 1.0);

    const bool dump = option.corrDump.length() > 0 && option.corrDump != "noDump";

    while (itr++ < option.maxItr) {
        for (int l = 0; l < option.dimension; ++l) {
            randPlacer(l) = randPair.generate();
        }
        randPlacer(3) = randPair_g.generate();
        auto tuple = state->gibbsSampling(randPlacer(0), randPlacer(1), randPlacer(2), randPlacer(3));
        double newHamiltonian = system.hamiltonian(tuple);

        double newEnergy = newHamiltonian, oldEnergy = hamiltonian;

        if (option.outputInt > 0 && itr % option.outputInt == 0 && itr / option.outputInt > option.skip) {
            cerr << oldEnergy << ' ' << state->getSum() << endl;
            if (dump) {
                state->updateCorr(corrMat_ii, corrMat_hh);
            }
        }
        if (coutInt > 0 && itr % (10*coutInt) == 0) {
            cout << format("\rProgress: % 3d%%; system energy: %.2f.") % (itr/coutInt) % oldEnergy;
        }

        if (newEnergy < oldEnergy) {
            state->flip(randPlacer(0), randPlacer(1), randPlacer(2), randPlacer(3), tuple);
            hamiltonian = newHamiltonian;
        } else {
            double alpha = unif_dis(unif_gen);
            double acceptProb = exp(-option.beta * (newEnergy - oldEnergy));
            if (alpha < acceptProb) {
                state->flip(randPlacer(0), randPlacer(1), randPlacer(2), randPlacer(3), tuple);
                hamiltonian = newHamiltonian;
            }
        }
    }

    if (dump) {
        ofstream out(option.corrDump.c_str(), std::ios::out);
        if (!out) {
            fprintf(stderr, "Cannot open the corr dump!\n");
            cout << endl;
            delete state;
            exit(4);
        }
        out << corrMat_ii << endl << corrMat_hh << endl;;
        out.close();
    }

    if (option.confDump.length() > 0 && option.confDump != "noDump") {
        ofstream out2(option.confDump.c_str(), std::ios::out);
        if (!out2) {
            fprintf(stderr, "Cannot open the conf dump!\n");
            cout << endl;
            delete state;
            exit(4);
        }
        state->print(out2);
        out2.close();
    }

    cout << endl;
    delete state;
    return 0;
}
