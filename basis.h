#ifndef __BASIS_H
#define __BASIS_H

#include <memory>

namespace itensor {

template <class Tensor>
struct Basis;

template <class Tensor>
using BasisPtr = std::unique_ptr<Basis<Tensor>>;

template <class Tensor>
struct Basis
    {
    Basis() { }
    virtual ~Basis() { }

    Tensor virtual
    state(int s, int n,
          const Args& args = Global::args()) const = 0;

    Tensor virtual
    newstate(int s, int n,
             const Args& args = Global::args()) const
        {
        return state(s,n,args);
        }

    Tensor virtual
    proj(int s, int n,
         const Args& args = Global::args()) const = 0;

    virtual
    const char*
    statestr(int s, int n,
             const Args& args = Global::args()) const = 0;

    };

};

#endif //__BASIS_H
