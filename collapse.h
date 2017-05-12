#ifndef __COLLAPSE_H
#define __COLLAPSE_H

#include "itensor/mps/mps.h"
#include "itensor/iqtensor.h"
#include "basis.h"
#include "itensor/all_basic.h"
#include "itensor/util/autovector.h"

namespace itensor {

template <typename Tensor>
std::vector<int>
collapse(MPSt<Tensor>& psi,
         const BasisPtr<Tensor>& B,
         const Args& args = Global::args())
    {
    const auto N = psi.N();
    const auto d = psi.sites()(1).m();
    //const auto d = psi.model()(1).m();

    std::vector<int> state(N+1);
        
    //note: in case of IQTensors, we need to convert the MPS iteratively into ITensor
    //otherwise symmetry breaking measurements (e.g. X-projections) cannot be applied!
    ITensor Aj1, Aj2;
       // Print(div(psi.A(1)));

        
    psi.position(1);
    for(int j = 1; j <= N; ++j)
        {
        //project into Itensor
        if(j==1) Aj1 = toITensor(psi.A(j));
        else     Aj1 = Aj2;
            
            
        Aj2 = toITensor(psi.A(j+1));

        //auto prob = Vector(d);
        std::vector<double> prob(d);
        Real tot = 0;
            
        //measure probabilities
        for(int s = 1; s < d; ++s)
            {
            //calculate overlap with projector
            auto z = (dag(prime(Aj1,Site))*toITensor(B->proj(j,s,args))*Aj1).cplx();
            prob[s-1] = (z.real());
                
            if(z.real() > 1.00000001 || z.real() < 0.  )
            {
                Print(z);
                Print(prob[s-1]);
                Error("projetor probability > 1 or < 0");
            }

            tot += prob[s-1];

            }
        prob[d-1] = 1-tot;

        //get a random number in [0,1]
        const auto r = Global::random();
            
#ifdef DEBUG
        if(r < 0 || r > 1) Error("Bad result from RNG: r outside interval [0,1]");
#endif

        //pick local product state
        auto& st = state.at(j);
        st = 1;
        Real pdisc = prob[st-1];

        while(r > pdisc)
            {
            ++st;
            pdisc += prob[st-1];
            }
            
        //project into product state
        psi.Anc(j) = B->newstate(j,st,args);
        if(j < N)
            {
            Aj2 *= dag(toITensor(B->state(j,st,args)))*Aj1;
            Aj2 *=1./std::sqrt(prob[st-1]);
            }
            
        }

    return state;
    }

};


#endif //__COLLAPSE_H
