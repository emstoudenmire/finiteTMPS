#ifndef __ROTATEXZ_H
#define __ROTATEXZ_H

#include "../basis.h"
#include "itensor/all_basic.h"

namespace itensor {

template <class Tensor>
BasisPtr<Tensor>
rotateXZ(const SiteSet& sites);

template <class Tensor>
struct RotateXZ : public Basis<Tensor>
    {
    RotateXZ(const SiteSet& sites)
        :
        N(sites.N()),
        projUpX(N+1),
        projDnX(N+1),
        projUpZ(N+1),
        projDnZ(N+1),
        stateUpX(N+1),
        stateDnX(N+1),
        stateUpZ(N+1),
        stateDnZ(N+1)
        {
        for(int j = 1; j <= N; ++j)
            {
            //indices (some redundancies)
            auto s_new = sites(j);
            auto s  = sites.si(j);
            auto sP = prime(s);
            auto Up = s(1);
            auto UpP = sP(1);
            auto Dn = s(2);
            auto DnP = sP(2);
            
            //projectors
            //const Tensor Id = sites.op("Id",j);
            projUpZ.at(j) = 0.5*sites.op("Id",j)+sites.op("Sz",j);
            projDnZ.at(j) = 0.5*sites.op("Id",j)-sites.op("Sz",j);

            //projUpX.at(j) = sites.op("Px+",j);
            auto& Pxu = projUpX.at(j);
            Pxu = mixedIQTensor(s,sP);
            Pxu.set(Up,UpP,+0.5);
            Pxu.set(Up,DnP,+0.5);
            Pxu.set(Dn,UpP,+0.5);
            Pxu.set(Dn,DnP,+0.5);

            //projDnX.at(j) = sites.op("Px-",j);
            auto& Pxd = projDnX.at(j);
            Pxd = mixedIQTensor(s,sP);
            Pxd.set(Up,UpP,+0.5);
            Pxd.set(Up,DnP,-0.5);
            Pxd.set(Dn,UpP,-0.5);
            Pxd.set(Dn,DnP,+0.5);



            IQTensor upX,
            dnX;
            IQTensor upZ(sites(j)),
            dnZ(sites(j));
            upX = mixedIQTensor(sites(j));//s,sP);
            dnX = mixedIQTensor(sites(j));//s,sP);

            
            IQIndexVal sup(sites(j,"Up")),
            sdn(sites(j,"Dn"));
            upZ = setElt(sup);
            dnZ = setElt(sdn);

            upX.set(s_new(1),ISqrt2);
            upX.set(s_new(2),ISqrt2);
            dnX.set(s_new(1),ISqrt2);
            dnX.set(s_new(2),-ISqrt2);

            //state projections
            stateUpZ.at(j) = upZ;
            stateDnZ.at(j) = dnZ;
            stateUpX.at(j) = upX;
            stateDnX.at(j) = dnX;
            }
        }

    Tensor
    state(int s, int n,
          const Args& args = Global::args()) const
        {
        return (n == 1 ? stateUpX.at(s) : stateDnX.at(s));
        }

    Tensor
    newstate(int s, int n,
             const Args& args = Global::args()) const
        {
        return (n == 1 ? stateUpZ.at(s) : stateDnZ.at(s));
        }

    Tensor
    proj(int s, int n,
         const Args& args = Global::args()) const
        {
        return (n == 1 ? projUpX.at(s) : projDnX.at(s));
        }

    const char*
    statestr(int s, int n,
             const Args& args = Global::args()) const
        { 
        return (n==1 ? "+" : "-"); 
        }

    private:

    int N;
    std::vector<Tensor> projUpX,
                        projDnX,
                        projUpZ,
                        projDnZ,
                        stateUpX,
                        stateDnX,
                        stateUpZ,
                        stateDnZ;

    };

template <class Tensor>
BasisPtr<Tensor>
rotateXZ(const SiteSet& sites)
    {
    return BasisPtr<Tensor>(new RotateXZ<Tensor>(sites));
    }

}; 

#endif 
