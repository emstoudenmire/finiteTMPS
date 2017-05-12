#ifndef __HEIS_OPS_H
#define __HEIS_OPS_H

#include "itensor/mps/siteset.h"

namespace itensor {

struct HeisOps
    {
    HeisOps(const SiteSet& sites,
            int Nx, int Ny,
            const Args& args = Global::args())
        :
        sites_(sites),
        putsite_(sites.N()+1,true),
        Nx_(Nx),
        Ny_(Ny)
        { 
        const auto J = args.getReal("J",1);
        Jz_ = args.getReal("Jz",J);
        Jxy_ = args.getReal("Jxy",J);
        Jx_ = args.getReal("Jx",0);
        hx_ = args.getReal("hx",0.);
        hz_ = args.getReal("hz",0.);
        //std::cout << boost::format("Jz_=%f,Jxy_=%f,hx_=%f,hz=%f") % Jz_% Jxy_%hx_%hz_ << std::endl;
        }

    IQTensor
    operator()(int i1, int i2, const std::string& type)
        {
        //std::cout << boost::format("(%d,%d)") % i1 % i2 << std::endl;
        IQTensor hh = sites_.op("Sz",i1)*sites_.op("Sz",i2) * Jz_;
        if(Jxy_ != 0)
            {
            hh += sites_.op("S+",i1)*sites_.op("S-",i2) * (Jxy_/2.);
            hh += sites_.op("S-",i1)*sites_.op("S+",i2) * (Jxy_/2.);
            }
        if(Jx_ != 0)
            {
            hh += sites_.op("Sx",i1)*sites_.op("Sx",i2) * (Jx_);
            }

        if(putsite_.at(i1))
            {
//            std::cout << "Adding field at " << i1 << std::endl;
            if(hx_ != 0) hh += sites_.op("Sx",i1)*sites_.op("Id",i2) * (-hx_);
            if(hz_ != 0) hh += sites_.op("Sz",i1)*sites_.op("Id",i2) * (-hz_);
            putsite_[i1] = false;
            }

        if(putsite_.at(i2))
            {
//            std::cout << "Adding field at " << i2 << std::endl;
            if(hx_ != 0) hh += sites_.op("Id",i1)*sites_.op("Sx",i2) * (-hx_);
            if(hz_ != 0) hh += sites_.op("Id",i1)*sites_.op("Sz",i2) * (-hz_);
            putsite_[i2] = false;
            }

        return hh;
        }

    private:

    const SiteSet& sites_;
    std::vector<bool> putsite_;
    int Nx_,
        Ny_;
    Real Jz_,
         Jxy_,
         Jx_,
         hz_,
         hx_;

    };

};

#endif
