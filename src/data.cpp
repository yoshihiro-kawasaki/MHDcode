#include "data.hpp"

Data::Data(const InputParameters &input)
{
    nx1_   = input.nx1_;
    nx2_   = input.nx2_;
    nx3_   = input.nx3_;
    x1min_ = input.x1min_;
    x1max_ = input.x1max_;
    x2min_ = input.x2min_;
    x2max_ = input.x2max_;
    x3min_ = input.x3min_;
    x3max_ = input.x3max_;
    time_.tlim_ = input.tlim_;
    time_.delta_tout_ = input.delta_tout_;
    dir_name_ = input.dir_name_;
    time_.t_ = 0.0;

    mkdir(dir_name_.c_str(), 0777);

    if (input.recon_type_ == ReconstructionType::DnonerCellSheme) {
        ngh_ = 1;
    } else if (input.recon_type_ == ReconstructionType::MusclMinmodScheme) {
        ngh_ = 2;
    }

    // set grid indices
    nx1totc_ = nx1_ + 2*ngh_;
    nx1totb_ = nx1_ + 2*ngh_ + 1;
    is_      = ngh_;
    ie_      = is_ + nx1_ - 1;
    xdim_    = DimensionsOfProblem::One;

    if (nx2_ > 1) {
        // 2D
        nx2totc_ = nx2_ + 2*ngh_;
        nx2totb_ = nx2_ + 2*ngh_ + 1;
        js_      = ngh_;
        je_      = js_ + nx2_ - 1;
        xdim_    = DimensionsOfProblem::Two;
    } else {
        // 1D
        nx2totc_ = 1;
        nx2totb_ = 1;
        js_      = 0;
        je_      = 0;
    }

    if (nx3_ > 1) {
        // 3D
        nx3totc_ = nx3_ + 2*ngh_;
        nx3totb_ = nx3_ + 2*ngh_ + 1;
        ks_      = ngh_;
        ke_      = ks_ + nx3_ - 1;
        xdim_    = DimensionsOfProblem::Three;
    } else {
        //2D
        nx3totc_ = 1;
        nx3totb_ = 1;
        ks_      = 0;
        ke_      = 0;
    }

    // allocate arrays
    q_    = array::Allocate4dArray<double>(NVAR, nx3totc_, nx2totc_, nx1totc_);
    u_    = array::Allocate4dArray<double>(NVAR, nx3totc_, nx2totc_, nx1totc_);
    
    x1c_  = array::Allocate1dArray<double>(nx1totc_);
    x2c_  = array::Allocate1dArray<double>(nx1totc_);
    x3c_  = array::Allocate1dArray<double>(nx1totc_);

    x1b_  = array::Allocate1dArray<double>(nx1totb_);
    x2b_  = array::Allocate1dArray<double>(nx2totb_);
    x3b_  = array::Allocate1dArray<double>(nx3totb_);

    flx1_ = array::Allocate4dArray<double>(NVAR, nx3totb_, nx2totb_, nx1totb_);
    flx2_ = array::Allocate4dArray<double>(NVAR, nx3totb_, nx2totb_, nx1totb_);
    flx3_ = array::Allocate4dArray<double>(NVAR, nx3totb_, nx2totb_, nx1totb_);

    // create grid

    dx1_ = (x1max_ - x1min_) / (nx1_);
    for (int i = 0; i < nx1totb_; ++i) {
        x1b_[i] = x1min_ + dx1_ * (i - (ngh_));
    }
    for (int i = 0; i < nx1totc_; ++i) {
        x1c_[i] = 0.5 * (x1b_[i] + x1b_[i+1]);
    }

    dx2_ = 0.0;
    if (xdim_ == DimensionsOfProblem::Two) {
        dx2_ = (x2max_ - x2min_) / nx2_;
        for (int i = 0; i < nx2totb_; ++i) {
            x2b_[i] = x2min_ + dx2_ * (i - (ngh_));
        }
        for (int i = 0; i < nx2totc_; ++i) {
            x2c_[i] = 0.5 * (x2b_[i] + x2b_[i+1]);
        }
    }

    dx3_ = 0.0;
    if (xdim_ == DimensionsOfProblem::Three) {
        dx3_ = (x3max_ - x3min_) / nx3_;
        for (int i = 0; i < nx3totb_; ++i) {
            x3b_[i] = x3min_ + dx3_ * (i - (ngh_));
        }
        for (int i = 0; i < nx3totc_; ++i) {
            x3c_[i] = 0.5 * (x3b_[i] + x3b_[i+1]);
        }
    }

    OutputGrid();
}

Data::~Data()
{
    if (is_allocate_array_) {
        array::Delete4dArray<double>(q_);
        array::Delete4dArray<double>(u_);

        array::Delete1dArray<double>(x1c_);
        array::Delete1dArray<double>(x2c_);
        array::Delete1dArray<double>(x3c_);

        array::Delete1dArray<double>(x1b_);
        array::Delete1dArray<double>(x2b_);
        array::Delete1dArray<double>(x3b_);

        array::Delete4dArray<double>(flx1_);
        array::Delete4dArray<double>(flx2_);
        array::Delete4dArray<double>(flx3_);

        is_allocate_array_ = false;
    }
}


void Data::Output(const long int count)
{
    int is = is_, ie = ie_, js = js_, je = je_, ks = ks_, ke = ke_;

    std::string filename = dir_name_ + "/st" + std::to_string(count);
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, filename);

    file << std::scientific;
    file << time_.t_ << std::endl;

    // for (int i = is; i <= ie; ++i) {
    //     file << x1c_[i] << " ";
    // }
    // file << std::endl;

    // for (int j = js; j <= je; ++j) {
    //     file << x2c_[j] << " ";
    // }
    // file << std::endl;

    // for (int k = ks; k <= ke; ++k) {
    //     file << x3c_[k] << " ";
    // }
    // file << std::endl;

    for (int k = ks; k <= ke; ++k) {
        for (int j = js; j <= je; ++j) {
            for (int i = is; i <= ie; ++i) {
                file << q_[IRHO][k][j][i] << " "
                     << q_[IV1][k][j][i] << " "
                     << q_[IV2][k][j][i] << " "
                     << q_[IV3][k][j][i] << " "
                     << q_[IPR][k][j][i] << " "
                     << q_[IB1][k][j][i] << " "
                     << q_[IB2][k][j][i] << " "
                     << q_[IB3][k][j][i] 
                     << std::endl; 
            }
        }
    }

    file.close();
    return;
}


void Data::OutputGrid()
{
    int is = is_, ie = ie_, js = js_, je = je_, ks = ks_, ke = ke_;

    std::string filename = dir_name_ + "/grid";
    std::ofstream file(filename, std::ios::trunc | std::ios::out);
    FAILED_TO_OPEN(file, filename);

    file << std::scientific;
    for (int i = is; i <= ie; ++i) {
        file << x1c_[i] << " ";
    }
    file << std::endl;

    for (int j = js; j <= je; ++j) {
        file << x2c_[j] << " ";
    }
    file << std::endl;


    for (int k = ks; k <= ke; ++k) {
        file << x3c_[k] << " ";
    }
    file << std::endl;
}