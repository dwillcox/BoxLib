
#undef  BL_LANG_CC
#ifndef BL_LANG_FORT
#define BL_LANG_FORT
#endif

#include "REAL.H"
#include "CONSTANTS.H"
#include "FLUXREG_F.H"
#include "ArrayLim.H"

#if defined(BL_USE_FLOAT) || defined(BL_T3E) || defined(BL_CRAY)
#define  BL_USE_FLOAT_CONST
#define  onetominuseight 1.0e-8
#else
#define  BL_USE_DOUBLE_CONST
#define  onetominuseight 1.0d-8
#endif

#define SDIM 3

c ::: -----------------------------------------------------------
c ::: Add fine grid flux to flux register.  Flux array is a fine grid
c ::: edge based object, Register is a coarse grid edge based object.      
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: DIMS(reg)  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: DIMS(flx)  => index limits for flx
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(3)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_FRFINEADD(reg,DIMS(reg),flx,DIMS(flx),
     &                          numcomp,dir,ratio,mult)

#ifdef BL_USE_OMP
      use omp_lib
#endif
      implicit none
      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    ratio(3), dir, numcomp
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)
      
      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

#ifdef BL_USE_OMP
      logical nested
      nested = omp_in_parallel()
#endif
      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)

      if (dir .eq. 0) then
c
c        ::::: flux normal to X direction
c
         ic = ARG_L1(reg)
         i = ic*ratiox
         if (ARG_L1(reg) .ne. ARG_H1(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (i .lt. ARG_L1(flx) .or. i .gt. ARG_H1(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

!$omp parallel if (.not. nested) private(n,koff,kc,k,joff,jc,j)
         do koff = 0, ratioz-1
!$omp do collapse(2)
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do joff = 0, ratioy-1            
                     do jc = ARG_L2(reg), ARG_H2(reg)
                        j = ratioy*jc + joff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
!$omp end do
         end do
!$omp end parallel

      else if (dir .eq. 1) then
c        ::::: flux normal to Y direction
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

!$omp parallel if (.not. nested) private(n,koff,kc,k,ioff,ic,i)
         do koff = 0, ratioz-1
!$omp do collapse(2)
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
!$omp end do
         end do
!$omp end parallel

      else
c
c        ::::: flux normal to Z direction
c
         kc = ARG_L3(reg)
         k = kc*ratioz
         if (ARG_L3(reg) .ne. ARG_H3(reg)) then
            call bl_abort("FORT_FRFINEADD: bad register direction")
         end if
         if (k .lt. ARG_L3(flx) .or. k .gt. ARG_H3(flx)) then
            call bl_abort("FORT_FRFINEADD: index outside flux range")
         end if

!$omp parallel if (.not. nested) private(n,joff,jc,j,ioff,ic,i)
         do joff = 0, ratioy-1
!$omp do collapse(2)
            do n = 1, numcomp
               do jc = ARG_L2(reg), ARG_H2(reg)
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + mult*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
!$omp end do
         end do
!$omp end parallel

      end if
      
      end

c ::: -----------------------------------------------------------
c ::: Add fine grid flux times area to flux register.  
c ::: Flux array is a fine grid edge based object, Register is a 
c ::: coarse grid edge based object.      
c ::: It is assumed that the coarsened flux region contains the register
c ::: region.
c ::: 
c ::: INPUTS/OUTPUTS:
c ::: reg       <=> edge centered coarse grid flux register
c ::: DIMS(reg)  => index limits for reg
c ::: flx        => edge centered fine grid flux array
c ::: DIMS(flx)  => index limits for flx
c ::: area       => edge centered area array
c ::: DIMS(area) => index limits for area
c ::: numcomp    => number of components to update
c ::: dir        => direction normal to flux register
c ::: ratio(3)   => refinement ratios between coarse and fine
c ::: mult       => scalar multiplicative factor      
c ::: -----------------------------------------------------------

      subroutine FORT_FRFAADD(reg,DIMS(reg),flx,DIMS(flx),area,DIMS(area),
     &                        numcomp,dir,ratio,mult)

#ifdef BL_USE_OMP
      use omp_lib
#endif
      implicit none
      integer    DIMDEC(reg)
      integer    DIMDEC(flx)
      integer    DIMDEC(area)
      integer    ratio(3), dir, numcomp
      REAL_T     mult
      REAL_T     reg(DIMV(reg),numcomp)
      REAL_T     flx(DIMV(flx),numcomp)
      REAL_T     area(DIMV(area))
      
      integer    n, i, j, k, ic, jc, kc, ioff, joff, koff
      integer    ratiox, ratioy, ratioz

#ifdef BL_USE_OMP
      logical nested
      nested = omp_in_parallel()
#endif
      ratiox = ratio(1)
      ratioy = ratio(2)
      ratioz = ratio(3)

      if (dir .eq. 0) then
c
c        ::::: flux normal to X direction
c
         ic = ARG_L1(reg)
         i = ic*ratiox
         if (ARG_L1(reg) .ne. ARG_H1(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (i .lt. ARG_L1(flx) .or. i .gt. ARG_H1(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

!$omp parallel if (.not. nested) private(n,koff,kc,k,joff,jc,j)
         do koff = 0, ratioz-1
!$omp do collapse(2)
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do joff = 0, ratioy-1            
                     do jc = ARG_L2(reg), ARG_H2(reg)
                        j = ratioy*jc + joff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + 
     &                       mult*area(i,j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
!$omp end do
         end do
!$omp end parallel

      else if (dir .eq. 1) then
c
c        ::::: flux normal to Y direction
c
         jc = ARG_L2(reg)
         j = jc*ratioy
         if (ARG_L2(reg) .ne. ARG_H2(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (j .lt. ARG_L2(flx) .or. j .gt. ARG_H2(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

!$omp parallel if (.not. nested) private(n,koff,kc,k,ioff,ic,i)
         do koff = 0, ratioz-1
!$omp do collapse(2)
            do n = 1, numcomp
               do kc = ARG_L3(reg), ARG_H3(reg)
                  k = ratioz*kc + koff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + 
     &                       mult*area(i,j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
!$omp end do
         end do
!$omp end parallel

      else
c
c        ::::: flux normal to Z direction
c
         kc = ARG_L3(reg)
         k = kc*ratioz
         if (ARG_L3(reg) .ne. ARG_H3(reg)) then
            call bl_abort("FORT_FRFAADD: bad register direction")
         end if
         if (k .lt. ARG_L3(flx) .or. k .gt. ARG_H3(flx)) then
            call bl_abort("FORT_FRFAADD: index outside flux range")
         end if

!$omp parallel if (.not. nested) private(n,joff,jc,j,ioff,ic,i)
         do joff = 0, ratioy-1
!$omp do collapse(2)
            do n = 1, numcomp
               do jc = ARG_L2(reg), ARG_H2(reg)
                  j = ratioy*jc + joff
                  do ioff = 0, ratiox-1            
                     do ic = ARG_L1(reg), ARG_H1(reg)
                        i = ratiox*ic + ioff
                        reg(ic,jc,kc,n) = reg(ic,jc,kc,n) + 
     &                       mult*area(i,j,k)*flx(i,j,k,n)
                     end do
                  end do
               end do
            end do
!$omp end do
         end do
!$omp end parallel

      end if
      
      end

c :: 
c :: --------------------------------------------------------------
c :: reflux:   reflux the data on the outer boundary of
c ::           a fine grid.  
c :: 
c :: Inputs/Outputs
c :: s           <=>  state data array
c :: DIMS(s)      =>  index limits of s array
c :: vol          =>  volume array 
c :: DIMS(vol)    =>  index limits of vol array
c :: reg          =>  flux register
c :: DIMS(reg)    =>  index limits of reg array
c :: lo,hi        =>  subregion of s array to be updated
c :: numcomp      =>  number of components to update
c :: mult         =>  multiplative factor (+1 or -1 depending on nomal)
c :: --------------------------------------------------------------
c :: 
      subroutine FORT_FRREFLUX (s,DIMS(s),vol,DIMS(vol),reg,DIMS(reg),
     &                          lo,hi,shft,numcomp,mult)

      implicit none

      integer    DIMDEC(s)
      integer    DIMDEC(vol)
      integer    DIMDEC(reg)
      integer    lo(SDIM), hi(SDIM), shft(SDIM)
      integer    numcomp
      REAL_T     mult
      REAL_T       s(DIMV(s),numcomp)
      REAL_T     vol(DIMV(vol))
      REAL_T     reg(DIMV(reg),numcomp)
      
      integer    n, i, j, k, ilo, jlo, klo

      ilo = shft(1)
      jlo = shft(2)
      klo = shft(3)

      do n = 1, numcomp
         do k = lo(3), hi(3)
            do j = lo(2), hi(2)
               do i = lo(1), hi(1)
                  s(i-ilo,j-jlo,k-klo,n) = s(i-ilo,j-jlo,k-klo,n)
     &                 + mult*reg(i,j,k,n)/vol(i-ilo,j-jlo,k-klo)
               end do
            end do
         end do
      end do
      
      end
