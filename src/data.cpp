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
    } else if (input.recon_type_ == ReconstructionType::MinmodScheme || 
               input.recon_type_ == ReconstructionType::VanLeerScheme) {
        ngh_ = 2;
    }

    // set grid indices
    nx1totv_ = nx1_ + 2*ngh_;
    nx1totf_ = nx1_ + 2*ngh_ + 1;
    is_      = ngh_;
    ie_      = is_ + nx1_ - 1;
    xdim_    = DimensionsOfProblem::One;

    if (nx2_ > 1) {
        // 2D
        nx2totv_ = nx2_ + 2*ngh_;
        nx2totf_ = nx2_ + 2*ngh_ + 1;
        js_      = ngh_;
        je_      = js_ + nx2_ - 1;
        xdim_    = DimensionsOfProblem::Two;
    } else {
        // 1D
        nx2totv_ = 1;
        nx2totf_ = 1;
        js_      = 0;
        je_      = 0;
    }

    if (nx3_ > 1) {
        // 3D
        nx3totv_ = nx3_ + 2*ngh_;
        nx3totf_ = nx3_ + 2*ngh_ + 1;
        ks_      = ngh_;
        ke_      = ks_ + nx3_ - 1;
        xdim_    = DimensionsOfProblem::Three;
    } else {
        //2D
        nx3totv_ = 1;
        nx3totf_ = 1;
        ks_      = 0;
        ke_      = 0;
    }

    // allocate arrays
    q_    = array::Allocate4dArray<double>(NVAR, nx3totv_, nx2totv_, nx1totv_);
    u_    = array::Allocate4dArray<double>(NVAR, nx3totv_, nx2totv_, nx1totv_);
    
    x1v_  = array::Allocate1dArray<double>(nx1totv_);
    x2v_  = array::Allocate1dArray<double>(nx2totv_);
    x3v_  = array::Allocate1dArray<double>(nx3totv_);
    dx1v_ = array::Allocate1dArray<double>(nx1totv_);
    dx2v_ = array::Allocate1dArray<double>(nx2totv_);
    dx3v_ = array::Allocate1dArray<double>(nx3totv_);

    x1f_  = array::Allocate1dArray<double>(nx1totf_);
    x2f_  = array::Allocate1dArray<double>(nx2totf_);
    x3f_  = array::Allocate1dArray<double>(nx3totf_);
    dx1f_ = array::Allocate1dArray<double>(nx1totv_);
    dx2f_ = array::Allocate1dArray<double>(nx2totv_);
    dx3f_ = array::Allocate1dArray<double>(nx3totv_);

    // create grid

    // x1-direction 
    dx1_ = (x1max_ - x1min_) / (nx1_);
    for (int i = 0; i < nx1totf_; ++i) {
        x1f_[i] = x1min_ + dx1_ * (i - (ngh_));
    }
    for (int i = 0; i < nx1totv_; ++i) {
        // dx1f_[i] = dx1_;
        dx1f_[i] = x1f_[i+1] - x1f_[i];
    }
    for (int i = 0; i < nx1totv_; ++i) {
        x1v_[i] = 0.5 * (x1f_[i] + x1f_[i+1]);
    }
    for (int i = 0; i < nx1totv_-1; ++i) {
        dx1v_[i] = x1v_[i+1] - x1v_[i-1];
    }

    dx2_ = 0.0;
    if (xdim_ == DimensionsOfProblem::Two) {
        dx2_ = (x2max_ - x2min_) / nx2_;
        for (int j = 0; j < nx2totf_; ++j) {
            x2f_[j] = x2min_ + dx2_ * (j - (ngh_));
        }
        for (int j = 0; j < nx1totv_; ++j) {
            // dx1f_[j] = dx2_;
            dx2f_[j] = x2f_[j+1] - x2f_[j];
        }
        for (int j = 0; j < nx2totv_; ++j) {
            x2v_[j] = 0.5 * (x2f_[j] + x2f_[j+1]);
        }
        for (int j = 0; j < nx2totv_-1; ++j) {
            dx2v_[j] = x2v_[j+1] - x2v_[j];
        }
    }

    dx3_ = 0.0;
    if (xdim_ == DimensionsOfProblem::Three) {
        dx3_ = (x3max_ - x3min_) / nx3_;
        for (int k = 0; k < nx3totf_; ++k) {
            x3f_[k] = x3min_ + dx3_ * (k - (ngh_));
        }
        for (int k = 0; k < nx3totv_; ++k) {
            dx3f_[k] = x3f_[k+1] - x3f_[k];
        }
        for (int k = 0; k < nx3totv_; ++k) {
            x3v_[k] = 0.5 * (x3f_[k] + x3f_[k+1]);
        }
        for (int k = 0; k < nx3totv_-1; ++k) {
            dx3v_[k] = x3v_[k+1] - x3v_[k];
        }
    }

    OutputGrid();
}

Data::~Data()
{
    if (is_allocate_array_) {
        array::Delete4dArray<double>(q_);
        array::Delete4dArray<double>(u_);

        array::Delete1dArray<double>(x1v_);
        array::Delete1dArray<double>(x2v_);
        array::Delete1dArray<double>(x3v_);
        array::Delete1dArray<double>(dx1v_);
        array::Delete1dArray<double>(dx2v_);
        array::Delete1dArray<double>(dx3v_);

        array::Delete1dArray<double>(x1f_);
        array::Delete1dArray<double>(x2f_);
        array::Delete1dArray<double>(x3f_);
        array::Delete1dArray<double>(dx1f_);
        array::Delete1dArray<double>(dx2f_);
        array::Delete1dArray<double>(dx3f_);

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
        file << x1v_[i] << " ";
    }
    file << std::endl;

    for (int j = js; j <= je; ++j) {
        file << x2v_[j] << " ";
    }
    file << std::endl;


    for (int k = ks; k <= ke; ++k) {
        file << x3v_[k] << " ";
    }
    file << std::endl;
}