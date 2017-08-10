#include "itensor/all.h"
#include "TStateObserver.h"
#include "S2.h"

using namespace std;
using namespace itensor;

using TensorT = IQTensor;
using MPOT = MPOt<TensorT>;
using MPST = MPSt<TensorT>;

int 
main(int argc, char* argv[])
    {
    printfln("TensorT == %s",(std::is_same<TensorT,ITensor>::value ? "ITensor" : "IQTensor"));

    //Get parameter file
    if(argc != 2)
        {
        printfln("Usage: %s inputfile.",argv[0]);
        return 0;
        }
    auto input = InputGroup(argv[1],"input");

    auto Nx = input.getInt("Nx",10);
    auto Ny = input.getInt("Ny",1);
    auto periodic = input.getYesNo("periodic",false);
    auto lattice_type = input.getString("lattice_type","triangular");

    auto beta = input.getReal("beta",1);
    auto tau = input.getReal("tau",0.05);

    auto maxm = input.getInt("maxm",1000);
    auto cutoff = input.getReal("cutoff",1E-11);

    auto Jz = input.getReal("Jz",1.);
    auto Jxy = input.getReal("Jxy",1.);

    auto realstep = input.getYesNo("realstep",false);
    auto verbose = input.getYesNo("verbose",false);

    auto N = Nx*Ny;

    Args args;
    args.add("Ny",Ny);
    args.add("Jz",Jz);
    args.add("Jxy",Jxy);
    args.add("Maxm",maxm);
    args.add("Cutoff",cutoff);
    args.add("YPeriodic",periodic);
    args.add("Verbose",verbose);

    auto sites = SpinHalf(2*N);
    writeToFile("sites",sites);

    LatticeGraph lattice; 
    if(lattice_type == "triangular")
        lattice = triangularLattice(Nx,Ny,args);
    else if(lattice_type == "square")
        lattice = squareLattice(Nx,Ny,args);

    auto ampo = AutoMPO(sites);
    for(auto b : lattice)
        {
        auto s1 = 2*b.s1-1,
             s2 = 2*b.s2-1;
        ampo += (0.5*Jxy),"S+",s1,"S-",s2;
        ampo += (0.5*Jxy),"S-",s1,"S+",s2;
        ampo +=        Jz,"Sz",s1,"Sz",s2;
        }

    MPOT expHa,expHb;
    MPOT expH;

    if(realstep)
        {
        expH = toExpH<TensorT>(ampo,tau);
        }
    else
        {
        auto taua = tau/2.*(1.+1._i);
        auto taub = tau/2.*(1.-1._i);
        println("Making expHa and expHb");
        expHa = toExpH<TensorT>(ampo,taua);
        expHb = toExpH<TensorT>(ampo,taub);
        }

    auto H = MPOT(ampo);

    auto S2 = makeS2(sites,{"SkipAncilla=",true});

    //
    // Make initial 'wavefunction' which is a product
    // of perfect singlets between neighboring sites
    //
    auto psi = MPST(sites);
    for(int n = 1; n <= 2*N; n += 2)
        {
        auto s1 = sites(n);
        auto s2 = sites(n+1);
        auto wf = TensorT(s1,s2);
        wf.set(s1(1),s2(2), ISqrt2);
        wf.set(s1(2),s2(1), -ISqrt2);
        TensorT D;
        psi.Aref(n) = TensorT(s1);
        psi.Aref(n+1) = TensorT(s2);
        svd(wf,psi.Aref(n),D,psi.Aref(n+1));
        psi.Aref(n) *= D;
        }

    auto obs = TStateObserver<TensorT>(psi);

    auto ttotal = beta/2.;
    const int nt = int(ttotal/tau+(1e-9*(ttotal/tau)));
    if(fabs(nt*tau-ttotal) > 1E-9)
        {
        Error("Timestep not commensurate with total time");
        }
    printfln("Doing %d steps of tau=%f",nt,tau);

    auto targs = args;

    auto En = Vector(nt);
    auto Sus = Vector(nt);
    auto Betas = Vector(nt);

    Real tsofar = 0;
    for(int tt = 1; tt <= nt; ++tt)
        {
        if(realstep)
            {
            psi = exactApplyMPO(expH,psi,args);
            }
        else
            {
            psi = exactApplyMPO(expHa,psi,args);
            psi = exactApplyMPO(expHb,psi,args);
            }
        psi.Aref(1) /= norm(psi.A(1));
        tsofar += tau;
        targs.add("TimeStepNum",tt);
        targs.add("Time",tsofar);
        targs.add("TotalTime",ttotal);
        obs.measure(targs);

        //Record beta value
        auto bb = (2*tsofar);
        Betas(tt-1) = bb;

        //
        // Measure Energy
        //
        auto en = overlap(psi,H,psi);
        printfln("\nEnergy/N %.4f %.20f",bb,en/N);
        En(tt-1) = en/N;

        //
        // Measure Susceptibility
        //
        auto s2val = overlap(psi,S2,psi);
        Sus(tt-1) = (s2val*bb/3.)/N;

        println();
        }

    std::ofstream enf("en.dat");
    std::ofstream susf("sus.dat");
    for(auto n : range(Betas))
        {
        enf << format("%.14f %.14f\n",Betas(n),En(n));
        susf << format("%.14f %.14f\n",Betas(n),Sus(n));
        }
    enf.close();
    susf.close();

    writeToFile("sites",sites);
    writeToFile("psi",psi);

    return 0;
    }

