#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>
#include <MultiFabUtil.H>

#include "Particles.H"

// declare routines below
void solve_for_accel(PArray<MultiFab>& rhs, PArray<MultiFab>& phi, PArray<MultiFab>& grad_phi, 
                     const Array<Geometry>& geom, int base_level, int finest_level, Real offset);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    // Define this as a 2-level problem.
    int nlevs = 2;

    // ********************************************************************************************
    // All of this defines the level 0 information -- size of box, type of boundary condition, etc.
    // ********************************************************************************************

    // This defines the physical size of the box.  Right now the box is [0,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       real_box.setLo(n,0.0);
       real_box.setHi(n,1.0);
    }

    // This defines the physical size of the fine box which we will use to constrain where we put the particles.
    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.4);
       fine_box.setHi(n,0.6);
    }

    // Define the lower and upper corner of a 3D domain
    IntVect domain_lo(0 , 0, 0); 
    int n_cell = 8;
    IntVect domain_hi(n_cell-1,n_cell-1,n_cell-1); 
 
    // Build a box for the level 0 domain
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Array<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // Now we make the refined level be the center eighth of the domain
    int n_fine = n_cell*rr[0];
    IntVect refined_lo(n_fine/4,n_fine/4,n_fine/4); 
    IntVect refined_hi(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1);

    // This says we are using Cartesian coordinates
    int coord = 0;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 

    // This defines a Geometry object which is useful for writing the plotfiles  
    Array<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, coord, is_per);
    for (int lev = 1; lev < nlevs; lev++)
    {
	geom[lev].define(BoxLib::refine(geom[lev-1].Domain(), rr[lev-1]),
			 &real_box, coord, is_per);
    }

    // ********************************************************************************************
    // This now defines the refinement information 
    // ********************************************************************************************

    // Build a box for the level 1 domain
    Box refined_patch(refined_lo, refined_hi);

    // Build an array of BoxArrays,
    // then initialize the level 0 BoxArray with the domain.
    Array<BoxArray> ba(nlevs);
    ba[0].define(domain);
    ba[1].define(refined_patch);

    // break the BoxArrays at both levels into 32^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(64);
    }

    // build a multifab for the rhs on the box array with 
    PArray<MultiFab> rhs; 
    PArray<MultiFab> phi;
    PArray<MultiFab> grad_phi;
    Array<DistributionMapping> dmap(nlevs);

    rhs.resize(nlevs,PArrayManage);
    phi.resize(nlevs,PArrayManage);
    grad_phi.resize(nlevs,PArrayManage);

    for (int lev = 0; lev < nlevs; lev++)
    {
	//                                    # componet  # ghost cells
	rhs.set     (lev,new MultiFab(ba[lev],1          ,0));
	phi.set     (lev,new MultiFab(ba[lev],1          ,1));
	grad_phi.set(lev,new MultiFab(ba[lev],BL_SPACEDIM,1));

	phi[lev].setVal(0.0);
	grad_phi[lev].setVal(0.0);

	dmap[lev] = rhs[lev].DistributionMap();
    }

    // Define a new particle container to hold my particles.
    // This holds a charge as well as three velocity components, three acceleration components  and three position components.
    typedef ParticleContainer<1+2*BL_SPACEDIM> MyParticleContainer;
    
    // Build a new particle container to hold my particles.
    MyParticleContainer* MyPC = new MyParticleContainer(geom,dmap,ba,rr);

    int num_particles = 1;
    int iseed = 10;
    Real mass  = 10.0;

    // Here we do a set of two experiments.  
    // 1) Do a single-level solve on level 0 as if there is no level 1, then 
    //    do a single-level solve on level 1 using the boundary conditions from level 0 from the previous step.
    // 2) Do a multi-level solve on levels 0 and 1.
    // We assume for all these tests that the particles within the area covered by the refined patch.

    // **************************************************************************
    // 1) Do a single-level solve on level 0, then do a solve on level 1 using the solution
    //    from level 0 for boundary condtions
    // **************************************************************************

    // Initialize "num_particles" number of particles, each with mass/charge "mass"
    bool serialize = false;
    MyPC->InitRandom(num_particles,iseed,mass,serialize,fine_box);
    MyPC->Redistribute();

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // **************************************************************************
    // Compute the total charge of all particles in order to compute the offset
    //     to make the Poisson equations solvable
    // **************************************************************************

    Real offset = 0.;
    if (geom[0].isAllPeriodic()) 
    {
        for (int lev = 0; lev < nlevs; lev++)
            offset = MyPC->sumParticleMass(lev);
        if (ParallelDescriptor::IOProcessor())
           std::cout << "Total charge of particles = " << offset << std::endl;
        offset /= geom[0].ProbSize();
    }

    // **************************************************************************

    // Define the density on level 0 from all particles at all levels
    int base_level   = 0;
    int finest_level = 1;

    PArray<MultiFab> PartMF;
    MyPC->AssignDensity(PartMF, base_level, 1, finest_level);

    for (int lev = finest_level - 1 - base_level; lev >= 0; lev--)
        BoxLib::average_down(PartMF[lev+1],PartMF[lev],0,1,rr[lev]);

    for (int lev = 0; lev < nlevs; lev++)
        MultiFab::Add(rhs[base_level+lev], PartMF[lev], 0, 0, 1, 0);
    // **************************************************************************
 
    // Define this to be solve at level 0 only
    base_level   = 0;
    finest_level = 0;

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solving for phi at level 0 ... " << std::endl;
    solve_for_accel(rhs,phi,grad_phi,geom,base_level,finest_level,offset);

    // Define this to be solve at level 1 only
    base_level   = 1;
    finest_level = 1;

    // Use multigrid to solve Lap(phi) = rhs with boundary conditions from level 0
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solving for phi at level 1 ... " << std::endl;
    solve_for_accel(rhs,phi,grad_phi,geom,base_level,finest_level,offset);
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Solved  for phi at level 1 ... " << std::endl;

    // Fill the particle data with the acceleration at the particle location
    // Note that we set dummy_dt = 0 so we don't actually move the particles.
    int start_comp = BL_SPACEDIM+1;
    Real dummy_dt = 0.0;
    for (int lev = 0; lev < nlevs; lev++)
        MyPC->moveKick(grad_phi[lev],0,dummy_dt,1.0,1.0,start_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after_level_solves");

    // **************************************************************************
    // 2) Do a multi-level solve on levels 0 and 1.
    // **************************************************************************

    // Reset everything to 0
    for (int lev = 0; lev < nlevs; lev++)
    {
	phi[lev].setVal(0.0);
	grad_phi[lev].setVal(0.0);
    }

    // Define this to be solve at multi-level solve
    base_level   = 0;
    finest_level = 1;

    // Redistribute the particles since we include both levels now.
    MyPC->Redistribute();

    // Use the PIC approach to deposit the "mass" onto the grid
    MyPC->AssignDensity(rhs,base_level,1,finest_level);

    // std::cout << "RHS " << rhs[0][0] << std::endl;

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_before");

    // Use multigrid to solve Lap(phi) = rhs with periodic boundary conditions (set above)
    solve_for_accel(rhs,phi,grad_phi,geom,base_level,finest_level,offset);

    // Fill the particle data with the acceleration at the particle location
    // Note that we set dummy_dt = 0 so we don't actually move the particles.
    start_comp = BL_SPACEDIM+1;
    dummy_dt = 0.0;
    MyPC->moveKick(grad_phi[0],0,dummy_dt,1.0,1.0,start_comp);

    // Write out the positions, masses and accelerations of each particle.
    MyPC->WriteAsciiFile("Particles_after_multilevel_solve");

    delete MyPC;

    BoxLib::Finalize();
}
