
#include <winstd.H>
#include <BArena.H>
#include <FluxRegister.H>
#include <Geometry.H>
#include <FLUXREG_F.H>
#include <ParallelDescriptor.H>
#include <BLProfiler.H>
#include <ccse-mpi.H>

#include <vector>

FluxRegister::FluxRegister ()
{
    fine_level = ncomp = -1;
    ratio = IntVect::TheUnitVector();
    ratio.scale(-1);
}

FluxRegister::FluxRegister (const BoxArray& fine_boxes, 
                            const IntVect&  ref_ratio,
                            int             fine_lev,
                            int             nvar)
{
    define(fine_boxes,ref_ratio,fine_lev,nvar);
}

FluxRegister::FluxRegister (const BoxArray&            fine_boxes, 
                            const IntVect&             ref_ratio,
                            int                        fine_lev,
                            int                        nvar,
                            const DistributionMapping& dm)
{
    define(fine_boxes,ref_ratio,fine_lev,nvar,dm);
}

const IntVect&
FluxRegister::refRatio () const
{
    return ratio;
}

int
FluxRegister::fineLevel () const
{
    return fine_level;
}

int
FluxRegister::crseLevel () const
{
    return fine_level-1;
}

int
FluxRegister::nComp () const
{
    return ncomp;
}

const BoxArray&
FluxRegister::coarsenedBoxes () const
{
    return grids;
}

void
FluxRegister::define (const BoxArray& fine_boxes, 
                      const IntVect&  ref_ratio,
                      int             fine_lev,
                      int             nvar)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

        IndexType typ(IndexType::TheCellType());

        typ.setType(dir,IndexType::NODE);

        BndryRegister::define(lo_face,typ,0,1,0,nvar);
        BndryRegister::define(hi_face,typ,0,1,0,nvar);
    }
}

void
FluxRegister::define (const BoxArray&            fine_boxes, 
                      const IntVect&             ref_ratio,
                      int                        fine_lev,
                      int                        nvar,
                      const DistributionMapping& dm)
{
    BL_ASSERT(fine_boxes.isDisjoint());
    BL_ASSERT(grids.size() == 0);

    ratio      = ref_ratio;
    fine_level = fine_lev;
    ncomp      = nvar;

    grids.define(fine_boxes);
    grids.coarsen(ratio);

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const Orientation lo_face(dir,Orientation::low);
        const Orientation hi_face(dir,Orientation::high);

        IndexType typ(IndexType::TheCellType());

        typ.setType(dir,IndexType::NODE);

        BndryRegister::define(lo_face,typ,0,1,0,nvar,dm);
        BndryRegister::define(hi_face,typ,0,1,0,nvar,dm);
    }
}

FluxRegister::~FluxRegister () {}

Real
FluxRegister::SumReg (int comp) const
{
    Real sum = 0.0;

    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        const FabSet& lofabs = bndry[Orientation(dir,Orientation::low) ];
        const FabSet& hifabs = bndry[Orientation(dir,Orientation::high)];

#ifdef _OPENMP
#pragma omp parallel reduction(+:sum)
#endif
        for (FabSetIter fsi(lofabs); fsi.isValid(); ++fsi)
        {
            sum += (lofabs[fsi].sum(comp) - hifabs[fsi].sum(comp));
        }
    }

    ParallelDescriptor::ReduceRealSum(sum);

    return sum;
}

struct Rec
{
    Rec (int         dIndex,
         int         sIndex,
         Orientation face,
         FillBoxId   fbid)
        :
        m_dIndex(dIndex),
        m_sIndex(sIndex),
        m_face(face),
        m_fbid(fbid) {}

    Rec (const IntVect& shift,
         int            dIndex,
         int            sIndex,
         Orientation    face,
         FillBoxId      fbid)
        :
        m_shift(shift),
        m_dIndex(dIndex),
        m_sIndex(sIndex),
        m_face(face),
        m_fbid(fbid) {}

    IntVect     m_shift;
    int         m_dIndex;
    int         m_sIndex;
    Orientation m_face;
    FillBoxId   m_fbid;
};

static
void
RefluxIt (const Rec&       rf,
          Real             scale,
          const Real*      multf,
          const BoxArray&  grids,
          MultiFab&        S,
          const MultiFab&  volume,
          const FArrayBox& reg,
          int              scomp,
          int              dcomp,
          int              ncomp)
{
    BL_ASSERT(S.DistributionMap()[rf.m_dIndex] == ParallelDescriptor::MyProc());
    BL_ASSERT(volume.DistributionMap()[rf.m_dIndex] == ParallelDescriptor::MyProc());

    Real mult;
    if (multf == 0)
        mult = rf.m_face.isLow() ? -scale : scale;
    else
        mult = (*multf)*scale;

    FArrayBox&       fab_S      = S[rf.m_dIndex];
    const FArrayBox& fab_volume = volume[rf.m_dIndex];
    Real*            s_dat      = fab_S.dataPtr(dcomp);
    const int*       slo        = fab_S.loVect();
    const int*       shi        = fab_S.hiVect();
    const Real*      vol_dat    = fab_volume.dataPtr();
    const Box&       fine_face  = BoxLib::adjCell(grids[rf.m_sIndex],rf.m_face);
    const Box&       sftbox     = S.box(rf.m_dIndex) + rf.m_shift;
    const Box&       ovlp       = sftbox & fine_face;
    const int*       lo         = ovlp.loVect();
    const int*       hi         = ovlp.hiVect();
    const int*       shft       = rf.m_shift.getVect();
    const int*       vlo        = fab_volume.loVect();
    const int*       vhi        = fab_volume.hiVect();
    const Real*      reg_dat    = reg.dataPtr(scomp);

    BL_ASSERT(ovlp.ok());

    FORT_FRREFLUX(s_dat,ARLIM(slo),ARLIM(shi),
                  vol_dat,ARLIM(vlo),ARLIM(vhi),
                  reg_dat,ARLIM(lo),ARLIM(hi),
                  lo,hi,shft,&ncomp,&mult);
}

void
FluxRegister::Reflux (MultiFab&       S,
                      const MultiFab& volume,
                      Real            scale,
                      int             src_comp,
                      int             dest_comp,
                      int             num_comp, 
                      const Geometry& geom,
		      const Real*     multf)
{
    BL_PROFILE("FluxRegister::Reflux()");

    FabSetId                          fsid[2*BL_SPACEDIM];
    std::vector<Rec>                  Recs;
    FabSetCopyDescriptor              fscd;
    std::vector< std::pair<int,Box> > isects;
    //
    // We use this to help "find" FluxRegisters with which we may intersect.
    // It assumes that FluxRegisters have width "1".
    //
    const int ng = 1;

    for (OrientationIter fi; fi; ++fi)
        fsid[fi()] = fscd.RegisterFabSet(&bndry[fi()]);

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int  idx = mfi.index();
        const Box& vbx = mfi.validbox();
        //
        // Find flux register that intersect with this grid.
        //
        grids.intersections(vbx,isects,ng);

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            const int k = isects[i].first;

            for (OrientationIter fi; fi; ++fi)
            {
                //
                // low (high) face of fine grid => high (low)
                // face of the exterior coarse grid cell updated.
                //
                const Orientation face = fi();

                const Box& ovlp = vbx & BoxLib::adjCell(grids[k],face);

                if (ovlp.ok())
                {
                    //
                    // Try to "minimize" the amount of data sent.
                    // Don't send the "whole" bndry register if not needed.
                    //
                    Box sbx = BoxLib::surroundingNodes(ovlp,face%BL_SPACEDIM);

                    sbx &= bndry[face].box(k);

                    FillBoxId fbid = fscd.AddBox(fsid[face],
                                                 sbx,
                                                 0,
                                                 k,
                                                 src_comp,
                                                 0,
                                                 num_comp);

                    Recs.push_back(Rec(idx,k,face,fbid));
                }
            }
        }
    }
    //
    // Add periodic possibilities.
    //
    if (geom.isAnyPeriodic())
    {
        Array<IntVect>  pshifts(27);

        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            const int  idx  = mfi.index();
            const Box& vbx  = mfi.validbox();

            for (int k = 0, N = grids.size(); k < N; k++)
            {
                const Box& bx = BoxLib::grow(grids[k],ng);

                if (geom.Domain().contains(bx)) continue;

                geom.periodicShift(bx,vbx,pshifts);

                const Box& kgrid = grids[k];

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv     = *it;
                    const Box&     sftbox = vbx + iv;

                    BL_ASSERT(bx.intersects(sftbox));

                    for (OrientationIter fi; fi; ++fi)
                    {
                        //
                        // low (high)  face of fine grid => high (low)
                        // face of the exterior coarse grid cell updated.
                        //
                        const Orientation face = fi();

                        const Box& ovlp = sftbox & BoxLib::adjCell(kgrid,face);

                        if (ovlp.ok())
                        {
                            //
                            // Try to "minimize" the amount of data sent.
                            // Don't send the "whole" bndry register if not needed.
                            //
                            Box sbx = BoxLib::surroundingNodes(ovlp,face%BL_SPACEDIM);

                            sbx &= bndry[face].box(k);

                            FillBoxId fbid = fscd.AddBox(fsid[face],
                                                         sbx,
                                                         0,
                                                         k,
                                                         src_comp,
                                                         0,
                                                         num_comp);

                            Recs.push_back(Rec(iv,idx,k,face,fbid));
                        }
                    }
                }
            }
        }
    }

    fscd.CollectData();

    const int N = Recs.size();

    for (int i = 0; i < N; i++)
    {
        const Rec&       rf   = Recs[i];
        const FillBoxId& fbid = rf.m_fbid;

        BL_ASSERT(bndry[rf.m_face].box(rf.m_sIndex).contains(fbid.box()));
        BL_ASSERT(S.DistributionMap()[rf.m_dIndex] == ParallelDescriptor::MyProc());
        BL_ASSERT(volume.DistributionMap()[rf.m_dIndex] == ParallelDescriptor::MyProc());

        FArrayBox reg(fbid.box(), num_comp);

        fscd.FillFab(fsid[rf.m_face], fbid, reg);

        RefluxIt(rf,scale,multf,grids,S,volume,reg,0,dest_comp,num_comp);
    }
}

void
FluxRegister::Reflux (MultiFab&       S,
                      Real            scale,
                      int             scomp,
                      int             dcomp,
                      int             ncomp, 
                      const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    MultiFab volume(S.boxArray(), 1, S.nGrow());

    volume.setVal(D_TERM(dx[0],*dx[1],*dx[2]), 0, 1, S.nGrow());

    Reflux(S,volume,scale,scomp,dcomp,ncomp,geom);
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        const MultiFab& area,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Orientation face_lo(dir,Orientation::low);
    const Orientation face_hi(dir,Orientation::high);
 
    MultiFab mf(mflx.boxArray(),numcomp,0);

#ifdef _OPENMP
#pragma omp parallel
#endif    
    for (MFIter mfi(mflx,true); mfi.isValid(); ++mfi)
    {
	const Box& bx = mfi.tilebox();
	
        mf[mfi].copy(mflx[mfi],bx,srccomp,bx,0,numcomp);

        mf[mfi].mult(mult,bx,0,numcomp);

        for (int i = 0; i < numcomp; i++)
            mf[mfi].mult(area[mfi],bx,bx,0,i,1);
    }

    for (int pass = 0; pass < 2; pass++)
    {
        const Orientation face = ((pass == 0) ? face_lo : face_hi);

        if (op == FluxRegister::COPY)
        {
            bndry[face].copyFrom(mf,0,0,destcomp,numcomp);
        }
        else
        {
            FabSet fs(bndry[face].boxArray(),numcomp);

            fs.setVal(0);

            fs.copyFrom(mf,0,0,0,numcomp);

            for (FabSetIter mfi(fs); mfi.isValid(); ++mfi)
                bndry[face][mfi].plus(fs[mfi],0,destcomp,numcomp);
        }
    }
}

void
FluxRegister::CrseInit (const MultiFab& mflx,
                        int             dir,
                        int             srccomp,
                        int             destcomp,
                        int             numcomp,
                        Real            mult,
                        FrOp            op)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= mflx.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    MultiFab area(mflx.boxArray(), 1, mflx.nGrow());

    area.setVal(1, 0, 1, area.nGrow());

    CrseInit(mflx,area,dir,srccomp,destcomp,numcomp,mult,op);
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],dir,k,srccomp,destcomp,numcomp,mult);
    }
}

void
FluxRegister::FineAdd (const MultiFab& mflx,
                       const MultiFab& area,
                       int             dir,
                       int             srccomp,
                       int             destcomp,
                       int             numcomp,
                       Real            mult)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mflx); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        FineAdd(mflx[mfi],area[k],dir,k,srccomp,destcomp,numcomp,mult);
    }
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Box&  flxbox = flux.box();
    const int*  flo    = flxbox.loVect();
    const int*  fhi    = flxbox.hiVect();
    const Real* flxdat = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

#ifndef NDEBUG
    Box cbox = BoxLib::coarsen(flux.box(),ratio);
    BL_ASSERT(cbox.contains(loreg.box()));
#endif
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFINEADD(lodat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

#ifndef NDEBUG
    BL_ASSERT(cbox.contains(hireg.box()));
#endif
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFINEADD(hidat,ARLIM(rlo),ARLIM(rhi),
                   flxdat,ARLIM(flo),ARLIM(fhi),
                   &numcomp,&dir,ratio.getVect(),&mult);
}

void
FluxRegister::FineAdd (const FArrayBox& flux,
                       const FArrayBox& area,
                       int              dir,
                       int              boxno,
                       int              srccomp,
                       int              destcomp,
                       int              numcomp,
                       Real             mult)
{
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= flux.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= ncomp);

    const Real* area_dat = area.dataPtr();
    const int*  alo      = area.loVect();
    const int*  ahi      = area.hiVect();
    const Box&  flxbox   = flux.box();
    const int*  flo      = flxbox.loVect();
    const int*  fhi      = flxbox.hiVect();
    const Real* flxdat   = flux.dataPtr(srccomp);

    FArrayBox& loreg = bndry[Orientation(dir,Orientation::low)][boxno];

#ifndef NDEBUG
    Box cbox = BoxLib::coarsen(flux.box(),ratio);
    BL_ASSERT(cbox.contains(loreg.box()));
#endif
    const int* rlo = loreg.box().loVect();
    const int* rhi = loreg.box().hiVect();
    Real* lodat = loreg.dataPtr(destcomp);
    FORT_FRFAADD(lodat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);

    FArrayBox& hireg = bndry[Orientation(dir,Orientation::high)][boxno];

#ifndef NDEBUG
    BL_ASSERT(cbox.contains(hireg.box()));
#endif
    rlo = hireg.box().loVect();
    rhi = hireg.box().hiVect();
    Real* hidat = hireg.dataPtr(destcomp);
    FORT_FRFAADD(hidat,ARLIM(rlo),ARLIM(rhi),
                 flxdat,ARLIM(flo),ARLIM(fhi),
                 area_dat,ARLIM(alo),ARLIM(ahi),
                 &numcomp,&dir,ratio.getVect(),&mult);
}

void
FluxRegister::write (const std::string& name, std::ostream& os) const
{
    if (ParallelDescriptor::IOProcessor())
    {
        os << ratio      << '\n';
        os << fine_level << '\n';
        os << ncomp      << '\n';
    }

    const BndryRegister* br = this;

    br->write(name,os);
}


void
FluxRegister::read (const std::string& name, std::istream& is)
{

    is >> ratio;
    is >> fine_level;
    is >> ncomp;

    BndryRegister* br = this;

    br->read(name,is);
}
