#ifndef __TROTTER_H
#define __TROTTER_H

#include <list>
#include "itensor/mps/bondgate.h"

namespace itensor {

template <typename Tensor>
using GateList = std::list<BondGate<Tensor>>;

template<class Tensor,class BondContainer,class OpFunction>
GateList<Tensor>
makeGates(const SiteSet& sites,
          const BondContainer& bonds,
          Real tau,
          OpFunction&& opf,
          const Args& args = Global::args());

template <typename GateT>
void
cleanGates(std::list<GateT>& gates);


//
// Implementations
//

template<class Tensor,class BondContainer,class OpFunction>
GateList<Tensor>
makeGates(const SiteSet& sites,
          const BondContainer& bonds,
          Real tau,
          OpFunction&& opf,
          const Args& args)
    {
    using GateT = BondGate<Tensor>;

    GateList<Tensor> gates;

    for(const auto& b : bonds)
        {
        const int i1 = b.s1,
                  i2 = b.s2; 

        Tensor hh = opf(i1,i2,b.type);

       // if(!hh.valid() || hh.norm() < 1E-12) continue;
        if(norm(hh) < 1E-12) continue;


        if(abs(i2-i1) == 1)
            {
            gates.push_back(GateT(sites,i1,i2,GateT::tImag,tau/2.,hh));

            }
        else
            {
            //Swap gates
            for(int k1 = i1; k1 <= i2-2; ++k1) 
                {

                gates.push_back(GateT(sites,k1,k1+1));
                }

            //Tensor that replaces index at i1 with index at i2-1
            Tensor II(sites.si(i1),dag(sites.si(i2-1)));
            for(int n = 1; n <= sites.si(i1).m(); ++n)
                {
                II.set(sites.si(i1)(n),sites.si(i2-1)(n),1.0);
                }

            hh *= II;
            hh *= dag(prime(II));

            gates.push_back(GateT(sites,i2-1,i2,GateT::tImag,tau/2.,hh));

            //Swap gates
            for(int k1 = i2-2; k1 >= i1; --k1) 
                {
                gates.push_back(GateT(sites,k1,k1+1));
                }
            }
        }

    // include b1.b2.b3....b3.b2.b1, second order trotter decomposition

    GateList<Tensor> gates2(gates);
    for(auto i = gates.rbegin(); i != gates.rend(); ++i)
        {
        gates2.push_back(*i);
        }
    gates = std::move(gates2);

    cleanGates(gates);

    return gates;
    }

template <typename GateT>
void
cleanGates(std::list<GateT>& gates)
    {
    bool success = true;
    while(success)
        {
        success = false;
        for(auto i = gates.begin(); i != gates.end(); ++i)
            {
            auto j = i;
            ++j;
            if(i->type() == GateT::Swap && j->type() == GateT::Swap && i->i1() == j->i1() && i->i2() == j->i2())
                {
                success = true;
                //std::cout << "erasing " << i->i() << " " << i->j() << std::endl;
                gates.erase(i,++j);
                break;
                }
            }
        }
    }

}; 

#endif //__TROTTER_H
