module cc_smoothers_module

  use bl_constants_module
  use cc_stencil_module
  use stencil_types_module

  implicit none

contains

  subroutine gs_line_solve_1d(ss, uu, ff, mm, lo, ng)

    use tridiag_module, only : tridiag

    integer, intent(in) :: lo(:)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in)    :: ff(lo(1):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:)
    real (kind = dp_t), intent(in)    :: ss(0:,lo(1):)
    integer            ,intent(in)    :: mm(lo(1):)

    real (kind = dp_t), allocatable :: a_ls(:), b_ls(:), c_ls(:), r_ls(:), u_ls(:)
    integer :: ilen, i, hi(size(lo))
    integer, parameter ::  XBC = 3

    hi = ubound(ff)

    ilen = hi(1)-lo(1)+1
    allocate(a_ls(0:hi(1)-lo(1)))
    allocate(b_ls(0:hi(1)-lo(1)))
    allocate(c_ls(0:hi(1)-lo(1)))
    allocate(r_ls(0:hi(1)-lo(1)))
    allocate(u_ls(0:hi(1)-lo(1)))

    do i = lo(1), hi(1)
      a_ls(i-lo(1)) = ss(2,i)
      b_ls(i-lo(1)) = ss(0,i)
      c_ls(i-lo(1)) = ss(1,i)
      r_ls(i-lo(1)) = ff(i)

      if ( hi(1) > lo(1) ) then
         if (bc_skewed(mm(i),1,+1)) then
            r_ls(i-lo(1)) = r_ls(i-lo(1)) - ss(XBC,i)*uu(i+2)
         else if (bc_skewed(mm(i),1,-1)) then
            r_ls(i-lo(1)) = r_ls(i-lo(1)) - ss(XBC,i)*uu(i-2)
         end if
      end if

    end do

    r_ls(0)           = r_ls(0)           - ss(2,lo(1)) * uu(lo(1)-1)
    r_ls(hi(1)-lo(1)) = r_ls(hi(1)-lo(1)) - ss(1,hi(1)) * uu(hi(1)+1)

    call tridiag(a_ls,b_ls,c_ls,r_ls,u_ls,ilen)
 
    do i = lo(1), hi(1)
       uu(i) = u_ls(i-lo(1))
    end do

  end subroutine gs_line_solve_1d

  subroutine gs_rb_smoother_1d(lo, hi, ss, uu, ff, mm, glo, ghi, ng, ns, n, skwd)
    integer,   intent(in) :: ng, ns, n
    integer,   intent(in) :: lo(1), hi(1), glo(1), ghi(1)
    real(dp_t),intent(in   )::ff(       glo(1)   :ghi(1))
    real(dp_t),intent(inout)::uu(       glo(1)-ng:ghi(1)+ng)
    real(dp_t),intent(in   )::ss(0:ns-1,glo(1)   :ghi(1))
    integer   ,intent(in   )::mm(       glo(1)   :ghi(1))
    logical,   intent(in) :: skwd

    integer :: i, looff, hioff
    real (kind = dp_t) :: dd(lo(1):hi(1))

    integer, parameter ::  XBC = 3
    real (kind = dp_t), parameter :: omega = 0.6_dp_t

    !! assumption: ss(0,i) many only vanish for 1x1 problems

    if ( glo(1) == ghi(1) ) then
       i = lo(1)
       if ( mod(i,2) == n ) then
          if ( ss(0,i) .eq. zero ) then
             dd(i) = ss(0,i)*uu(i) &
                  + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1) 
             uu(i) = uu(i) + (omega/ss(0,i)) * (ff(i) - dd(i))
          end if
       end if
       return
    end if

    if (mod(lo(1),2) .eq. n) then
       looff = 0
    else
       looff = 1
    end if
    if (mod(hi(1),2) .eq. n) then
       hioff = 0
    else
       hioff = 1
    end if

    do i = lo(1)+looff, hi(1)-hioff, 2
       dd(i) = ss(0,i)*uu(i) &
            + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1) 
    end do

    if (skwd) then
       i = lo(1)+looff
       if (i .eq. glo(1)) then
          if (bc_skewed(mm(i),1,+1)) then
             dd(i) = dd(i) + ss(XBC,i)*uu(i+2)
          end if
       end if

       i = hi(1)-hioff
       if (i .eq. ghi(1)) then
          if (bc_skewed(mm(i),1,-1)) then
             dd(i) = dd(i) + ss(XBC,i)*uu(i-2)
          end if
       end if
    end if

    do i = lo(1)+looff, hi(1)-hioff, 2
       uu(i) = uu(i) + (omega/ss(0,i)) * (ff(i) - dd(i))
    end do

  end subroutine gs_rb_smoother_1d

  subroutine gs_rb_smoother_2d(lo, hi, ss, uu, ff, mm, glo, ghi, ng, ns, n, skwd)
    integer,   intent(in) :: ng, ns, n
    integer,   intent(in) :: lo(2), hi(2), glo(2), ghi(2)
    real(dp_t),intent(in   )::ff(       glo(1)   :ghi(1)   ,glo(2)   :ghi(2)   )
    real(dp_t),intent(inout)::uu(       glo(1)-ng:ghi(1)+ng,glo(2)-ng:ghi(2)+ng)
    real(dp_t),intent(in   )::ss(0:ns-1,glo(1)   :ghi(1)   ,glo(2)   :ghi(2)   )
    integer   ,intent(in   )::mm(       glo(1)   :ghi(1)   ,glo(2)   :ghi(2)   )
    logical,   intent(in) :: skwd

    integer :: j, i, looff, hioff
    real (kind = dp_t) :: dd(lo(1):hi(1))

    integer, parameter ::  XBC = 5, YBC = 6

    !! assumption: ss(0,i,j) may only vanishe for 1x1 problems
    if ( all(glo == ghi) ) then
       i = lo(1); j = lo(2)
       if ( mod(i + j,2) == n ) then
          if ( ss(0,i,j) .eq. zero ) then
             dd(i) = ss(0,i,j)*uu(i,j) &
                  + ss(1,i,j)*uu(i+1,j) + ss(2,i,j)*uu(i-1,j) &
                  + ss(3,i,j)*uu(i,j+1) + ss(4,i,j)*uu(i,j-1)
             uu(i,j) = uu(i,j) + (one/ss(0,i,j))*(ff(i,j) - dd(i)) 
          end if
       end if
       return
    end if

    if (mod(lo(1)+lo(2),2) .eq. n) then
       looff = 0
    else
       looff = 1
    end if
    if (mod(hi(1)+lo(2),2) .eq. n) then
       hioff = 0
    else
       hioff = 1
    end if

    looff = 1-looff
    hioff = 1-hioff

    do j = lo(2),hi(2)
       looff = 1-looff
       hioff = 1-hioff

       do i = lo(1) + looff, hi(1) - hioff, 2
          dd(i) = ss(0,i,j)*uu(i,j) &
               + ss(1,i,j) * uu(i+1,j) + ss(2,i,j) * uu(i-1,j) &
               + ss(3,i,j) * uu(i,j+1) + ss(4,i,j) * uu(i,j-1)
       end do

       if (skwd) then
          i = lo(1) + looff
          if (i .eq. glo(1)) then
             if (bc_skewed(mm(i,j),1,+1)) then
                dd(i) = dd(i) + ss(XBC,i,j)*uu(i+2,j)
             end if
          end if

          i = hi(1) - hioff
          if (i .eq. ghi(1)) then
             if (bc_skewed(mm(i,j),1,-1)) then
                dd(i) = dd(i) + ss(XBC,i,j)*uu(i-2,j)
             end if
          end if

          if (j .eq. glo(2)) then
             do i = lo(1) + looff, hi(1) - hioff, 2
                if (bc_skewed(mm(i,j),2,+1)) then
                   dd(i) = dd(i) + ss(YBC,i,j)*uu(i,j+2)
                end if
             end do
          end if

          if (j .eq. ghi(2)) then
             do i = lo(1) + looff, hi(1) - hioff, 2
                if (bc_skewed(mm(i,j),2,-1)) then
                   dd(i) = dd(i) + ss(YBC,i,j)*uu(i,j-2)
                end if
             end do
          end if
       end if

       do i = lo(1) + looff, hi(1) - hioff, 2
          uu(i,j) = uu(i,j) + (one/ss(0,i,j))*(ff(i,j) - dd(i)) 
       end do
    end do

  end subroutine gs_rb_smoother_2d

  subroutine gs_rb_smoother_3d(lo, hi, ss, uu, ff, mm, glo, ghi, ng, ns, n, skwd)
    integer,   intent(in) :: ng, ns, n
    integer,   intent(in) :: lo(3), hi(3), glo(3), ghi(3)
    real(dp_t),intent(in   )::ff(       glo(1)   :ghi(1)   ,glo(2)   :ghi(2)   ,glo(3)   :ghi(3)   )
    real(dp_t),intent(inout)::uu(       glo(1)-ng:ghi(1)+ng,glo(2)-ng:ghi(2)+ng,glo(3)-ng:ghi(3)+ng)
    real(dp_t),intent(in   )::ss(0:ns-1,glo(1)   :ghi(1)   ,glo(2)   :ghi(2)   ,glo(3)   :ghi(3)   )
    integer   ,intent(in   )::mm(       glo(1)   :ghi(1)   ,glo(2)   :ghi(2)   ,glo(3)   :ghi(3)   )
    logical,   intent(in) :: skwd

    integer    :: i, j, k, looff, hioff
    real(dp_t) :: dd(lo(1):hi(1))

    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9
    real(dp_t), parameter :: omega = 1.15_dp_t

    if ( all(glo == ghi) ) then
       k = lo(3); j = lo(2); i = lo(1)
       if ( mod(i + j + k, 2) == n ) then
          if (abs(ss(0,i,j,k)) .gt. zero) then
             dd(i) = ss(0,i,j,k)*uu(i,j,k)   + &
                  ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                  ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                  ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)
             uu(i,j,k) = uu(i,j,k) + (omega/ss(0,i,j,k))*(ff(i,j,k) - dd(i))
          end if
       end if
       return
    end if

    do k = lo(3), hi(3)

       if (mod(lo(1)+lo(2)+k,2) .eq. n) then
          looff = 0
       else
          looff = 1
       end if
       if (mod(hi(1)+lo(2)+k,2) .eq. n) then
          hioff = 0
       else
          hioff = 1
       end if
       
       looff = 1-looff
       hioff = 1-hioff

       do j = lo(2), hi(2)
          looff = 1-looff
          hioff = 1-hioff

          do i = lo(1) + looff, hi(1) - hioff, 2
             dd(i) = ss(0,i,j,k)*uu(i,j,k)   + &
                  &  ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                  &  ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                  &  ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)
          end do
             
          if (skwd) then
             i = lo(1) + looff
             if (i .eq. glo(1)) then
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd(i) = dd(i) + ss(XBC,i,j,k)*uu(i+2,j,k)
                end if
             end if

             i = hi(1) - hioff
             if (i .eq. ghi(1)) then
                if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd(i) = dd(i) + ss(XBC,i,j,k)*uu(i-2,j,k)
                end if
             end if

             if (j .eq. glo(2)) then
                do i = lo(1) + looff, hi(1) - hioff, 2
                   if (bc_skewed(mm(i,j,k),2,+1)) then
                      dd(i) = dd(i) + ss(YBC,i,j,k)*uu(i,j+2,k)
                   end if
                end do
             end if

             if (j .eq. ghi(2)) then
                do i = lo(1) + looff, hi(1) - hioff, 2
                   if (bc_skewed(mm(i,j,k),2,-1)) then
                      dd(i) = dd(i) + ss(YBC,i,j,k)*uu(i,j-2,k)
                   end if
                end do
             end if

             if (k .eq. glo(3)) then
                do i = lo(1) + looff, hi(1) - hioff, 2
                   if (bc_skewed(mm(i,j,k),3,+1)) then
                      dd(i) = dd(i) + ss(ZBC,i,j,k)*uu(i,j,k+2)
                   end if
                end do
             end if

             if (k .eq. ghi(3)) then
                do i = lo(1) + looff, hi(1) - hioff, 2
                   if (bc_skewed(mm(i,j,k),3,-1)) then
                      dd(i) = dd(i) + ss(ZBC,i,j,k)*uu(i,j,k-2)
                   end if
                end do
             end if
          end if

          do i = lo(1) + looff, hi(1) - hioff, 2
             uu(i,j,k) = uu(i,j,k) + (omega/ss(0,i,j,k))*(ff(i,j,k) - dd(i))
          end do
       end do
    end do

  end subroutine gs_rb_smoother_3d

  subroutine fourth_order_smoother_2d(ss, uu, ff, lo, ng, stencil_type, n)
    use bl_prof_module
    integer           , intent(in) :: ng, n
    integer           , intent(in) :: lo(:)
    real (kind = dp_t), intent(in) :: ff(lo(1):, lo(2):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:,lo(1):, lo(2):)
    integer           , intent(in   ) :: stencil_type

    integer            :: i, j, hi(size(lo)), ioff
    real (kind = dp_t) :: dd

    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("fourth_order_smoother_2d")
    call build(bpt, "fourth_order_smoother_2d")

    hi = ubound(ff)

    if (stencil_type .eq. HO_CROSS_STENCIL) then

       do j = lo(2),hi(2)
          ioff = 0; if ( mod(lo(1) + j, 2) /= n ) ioff = 1
          do i = lo(1) + ioff, hi(1), 2
             if (abs(ss(0,i,j)) .gt. zero) then
               dd =   ss(0,i,j) * uu(i,j) &
                    + ss(1,i,j) * uu(i-2,j) + ss(2,i,j) * uu(i-1,j) &
                    + ss(3,i,j) * uu(i+1,j) + ss(4,i,j) * uu(i+2,j) &
                    + ss(5,i,j) * uu(i,j-2) + ss(6,i,j) * uu(i,j-1) &
                    + ss(7,i,j) * uu(i,j+1) + ss(8,i,j) * uu(i,j+2)
               uu(i,j) = uu(i,j) + (one/ss(0,i,j)) * (ff(i,j) - dd)
             end if
          end do
       end do

    else if (stencil_type .eq. HO_DENSE_STENCIL) then

       do j = lo(2),hi(2)
          do i = lo(1), hi(1)
             if ( abs(ss(0,i,j)) .gt. zero ) then
               dd =   ss( 0,i,j) * uu(i,j) &
                    + ss( 1,i,j) * uu(i-2,j-2) + ss( 2,i,j) * uu(i-1,j-2) & ! AT J-2
                    + ss( 3,i,j) * uu(i  ,j-2) + ss( 4,i,j) * uu(i+1,j-2) & ! AT J-2
                    + ss( 5,i,j) * uu(i+2,j-2)                            & ! AT J-2
                    + ss( 6,i,j) * uu(i-2,j-1) + ss( 7,i,j) * uu(i-1,j-1) & ! AT J-1
                    + ss( 8,i,j) * uu(i  ,j-1) + ss( 9,i,j) * uu(i+1,j-1) & ! AT J-1
                    + ss(10,i,j) * uu(i+2,j-1)                            & ! AT J-1
                    + ss(11,i,j) * uu(i-2,j  ) + ss(12,i,j) * uu(i-1,j  ) & ! AT J
                    + ss(13,i,j) * uu(i+1,j  ) + ss(14,i,j) * uu(i+2,j  ) & ! AT J
                    + ss(15,i,j) * uu(i-2,j+1) + ss(16,i,j) * uu(i-1,j+1) & ! AT J+1
                    + ss(17,i,j) * uu(i  ,j+1) + ss(18,i,j) * uu(i+1,j+1) & ! AT J+1
                    + ss(19,i,j) * uu(i+2,j+1)                            & ! AT J+1
                    + ss(20,i,j) * uu(i-2,j+2) + ss(21,i,j) * uu(i-1,j+2) & ! AT J+2
                    + ss(22,i,j) * uu(i  ,j+2) + ss(23,i,j) * uu(i+1,j+2) & ! AT J+2
                    + ss(24,i,j) * uu(i+2,j+2)                              ! AT J+2

               uu(i,j) = uu(i,j) + (one/ss(0,i,j)) * (ff(i,j) - dd)

             end if
          end do
       end do

    end if

    call destroy(bpt)
    call bl_proffortfuncstop("fourth_order_smoother_2d")

  end subroutine fourth_order_smoother_2d

  subroutine fourth_order_smoother_3d(ss, uu, ff, lo, ng, stencil_type)
    use bl_prof_module
    integer           , intent(in   ) :: ng
    integer           , intent(in   ) :: lo(:)
    real (kind = dp_t), intent(in   ) :: ff(lo(1):, lo(2):, lo(3):)
    real (kind = dp_t), intent(inout) :: uu(lo(1)-ng:, lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: ss(0:, lo(1):, lo(2):, lo(3):)
    integer           , intent(in   ) :: stencil_type

    integer            :: i, j, k, hi(size(lo))
    real (kind = dp_t) :: dd

    type(bl_prof_timer), save :: bpt

    call bl_proffortfuncstart("fourth_order_smoother_3d")
    call build(bpt, "fourth_order_smoother_3d")

    hi = ubound(ff)

    ! This is the fourth order stencil for constant coefficients.
    if (stencil_type .eq. HO_CROSS_STENCIL) then

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1), hi(1)
             if ( abs(ss(0,i,j,k)) .gt. zero ) then
               dd =   ss( 0,i,j,k) * uu(i,j,k) &
                    + ss( 1,i,j,k) * uu(i-2,j,k) + ss( 2,i,j,k) * uu(i-1,j,k) &
                    + ss( 3,i,j,k) * uu(i+1,j,k) + ss( 4,i,j,k) * uu(i+2,j,k) &
                    + ss( 5,i,j,k) * uu(i,j-2,k) + ss( 6,i,j,k) * uu(i,j-1,k) &
                    + ss( 7,i,j,k) * uu(i,j+1,k) + ss( 8,i,j,k) * uu(i,j+2,k) &
                    + ss( 9,i,j,k) * uu(i,j,k-2) + ss(10,i,j,k) * uu(i,j,k-1) &
                    + ss(11,i,j,k) * uu(i,j,k+1) + ss(12,i,j,k) * uu(i,j,k+2)
               uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd)
             end if
          end do
       end do
       end do

    ! This is the fourth order stencil for variable coefficients.
    else if (stencil_type .eq. HO_DENSE_STENCIL) then

       do k = lo(3),hi(3)
       do j = lo(2),hi(2)
          do i = lo(1), hi(1)
             if (abs(ss(0,i,j,k)) .gt. zero) then

                dd = &
                       ss( 0,i,j,k) * uu(i  ,j  ,k  ) &
                       ! Contributions from k-2
                     + ss( 1,i,j,k) * uu(i  ,j-2,k-2) + ss( 2,i,j,k) * uu(i  ,j-1,k-2) &
                     + ss( 3,i,j,k) * uu(i-2,j  ,k-2) + ss( 4,i,j,k) * uu(i-1,j  ,k-2) &
                     + ss( 5,i,j,k) * uu(i  ,j  ,k-2) + ss( 6,i,j,k) * uu(i+1,j  ,k-2) &
                     + ss( 7,i,j,k) * uu(i+2,j  ,k-2) + ss( 8,i,j,k) * uu(i  ,j+1,k-2) &
                     + ss( 9,i,j,k) * uu(i  ,j+2,k-2)                                  &
                       ! Contributions from k-1
                     + ss(10,i,j,k) * uu(i  ,j-2,k-1) + ss(11,i,j,k) * uu(i  ,j-1,k-1) &
                     + ss(12,i,j,k) * uu(i-2,j  ,k-1) + ss(13,i,j,k) * uu(i-1,j  ,k-1) &
                     + ss(14,i,j,k) * uu(i  ,j  ,k-1) + ss(15,i,j,k) * uu(i+1,j  ,k-1) &
                     + ss(16,i,j,k) * uu(i+2,j  ,k-1) + ss(17,i,j,k) * uu(i  ,j+1,k-1) &
                     + ss(18,i,j,k) * uu(i  ,j+2,k-1)                                  &
                       ! Contributions from j-2,k
                     + ss(19,i,j,k) * uu(i-2,j-2,k  ) + ss(20,i,j,k) * uu(i-1,j-2,k  ) &
                     + ss(21,i,j,k) * uu(i  ,j-2,k  ) + ss(22,i,j,k) * uu(i+1,j-2,k  ) &
                     + ss(23,i,j,k) * uu(i+2,j-2,k  )

                dd = dd &
                       ! Contributions from j-1,k
                     + ss(24,i,j,k) * uu(i-2,j-1,k  ) + ss(25,i,j,k) * uu(i-1,j-1,k  ) &
                     + ss(26,i,j,k) * uu(i  ,j-1,k  ) + ss(27,i,j,k) * uu(i+1,j-1,k  ) &
                     + ss(28,i,j,k) * uu(i+2,j-1,k  )                                  &
                       ! Contributions from j  ,k
                     + ss(29,i,j,k) * uu(i-2,j  ,k  ) + ss(30,i,j,k) * uu(i-1,j  ,k  ) &
                                                      + ss(31,i,j,k) * uu(i+1,j  ,k  ) &
                     + ss(32,i,j,k) * uu(i+2,j  ,k  )                                  &
                       ! Contributions from j+1,k
                     + ss(33,i,j,k) * uu(i-2,j+1,k  ) + ss(34,i,j,k) * uu(i-1,j+1,k  ) &
                     + ss(35,i,j,k) * uu(i  ,j+1,k  ) + ss(36,i,j,k) * uu(i+1,j+1,k  ) &
                     + ss(37,i,j,k) * uu(i+2,j+1,k  )                                  &
                       ! Contributions from j+2,k
                     + ss(38,i,j,k) * uu(i-2,j+2,k  ) + ss(39,i,j,k) * uu(i-1,j+2,k  ) &
                     + ss(40,i,j,k) * uu(i  ,j+2,k  ) + ss(41,i,j,k) * uu(i+1,j+2,k  ) &
                     + ss(42,i,j,k) * uu(i+2,j+2,k  )

                dd = dd &
                       ! Contributions from k+1
                     + ss(43,i,j,k) * uu(i  ,j-2,k+1) + ss(44,i,j,k) * uu(i  ,j-1,k+1) &
                     + ss(45,i,j,k) * uu(i-2,j  ,k+1) + ss(46,i,j,k) * uu(i-1,j  ,k+1) &
                     + ss(47,i,j,k) * uu(i  ,j  ,k+1) + ss(48,i,j,k) * uu(i+1,j  ,k+1) &
                     + ss(49,i,j,k) * uu(i+2,j  ,k+1) + ss(50,i,j,k) * uu(i  ,j+1,k+1) &
                     + ss(51,i,j,k) * uu(i  ,j+2,k+1)                                  &
                       ! Contributions from k+2
                     + ss(52,i,j,k) * uu(i  ,j-2,k+2) + ss(53,i,j,k) * uu(i  ,j-1,k+2) &
                     + ss(54,i,j,k) * uu(i-2,j  ,k+2) + ss(55,i,j,k) * uu(i-1,j  ,k+2) &
                     + ss(56,i,j,k) * uu(i  ,j  ,k+2) + ss(57,i,j,k) * uu(i+1,j  ,k+2) &
                     + ss(58,i,j,k) * uu(i+2,j  ,k+2) + ss(59,i,j,k) * uu(i  ,j+1,k+2) &
                     + ss(60,i,j,k) * uu(i  ,j+2,k+2)

               uu(i,j,k) = uu(i,j,k) + (one/ss(0,i,j,k)) * (ff(i,j,k) - dd)
             end if
          end do
       end do
       end do

    end if

    call destroy(bpt)
    call bl_proffortfuncstop("fourth_order_smoother_3d")

  end subroutine fourth_order_smoother_3d

  subroutine jac_smoother_1d(ss, uu, ff, mm, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  ss(0:,:)
    real (kind = dp_t), intent(in) ::  ff(:)
    integer           , intent(in) ::  mm(:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:)
    real (kind = dp_t):: wrk(size(ff,1))
    integer :: nx
    integer :: i
    real (kind = dp_t) :: dd
    integer, parameter ::  XBC = 3

    nx = size(ff,dim=1)
    do i = 1, nx
       dd = + ss(1,i)*uu(i+1) + ss(2,i)*uu(i-1) 
       if ( nx > 1 ) then
          if ( bc_skewed(mm(i),1,+1) ) then
             dd = dd + ss(XBC,i)*uu(i+2)
          else if ( bc_skewed(mm(i),1,-1) ) then
             dd = dd + ss(XBC,i)*uu(i-2)
          end if
       end if
       wrk(i) = ff(i) - dd
    end do
    do i = 1, nx
       uu(i) = uu(i) + (wrk(i)/ss(0,i)-uu(i))
    end do

  end subroutine jac_smoother_1d

  subroutine jac_dense_smoother_1d(ss, uu, ff, ng)
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) ::  ss(:,:)
    real (kind = dp_t), intent(in) ::  ff(:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:)
    real (kind = dp_t):: wrk(size(ff,1))
    integer :: nx, i

    nx = size(ff,dim=1)
    do i = 1, nx
       wrk(i) = ( ff(i) - ss(1,i)*uu(i+1) - ss(3,i)*uu(i-1) )
    end do
    do i = 1, nx
       uu(i) = uu(i) + (wrk(i)/ss(2,i)-uu(i))
    end do

  end subroutine jac_dense_smoother_1d

  subroutine jac_smoother_2d(ss, uu, ff, mm, ng)
    use bl_prof_module
    integer, intent(in) :: ng

    real (kind = dp_t), intent(in) ::  ss(0:,:,:)
    real (kind = dp_t), intent(in) ::  ff(:,:)
    integer           , intent(in) ::  mm(:,:)
    real (kind = dp_t), intent(inout) ::  uu(1-ng:,1-ng:)
    real (kind = dp_t) :: wrk(size(ff,1),size(ff,2))
    integer :: nx, ny
    integer :: i, j
    real (kind = dp_t) :: dd
    integer, parameter ::  XBC = 5, YBC = 6

    real (kind = dp_t), parameter :: omega = 4.0_dp_t/5.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_smoother_2d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)

    do j = 1, ny
       do i = 1, nx
          dd =    + ss(1,i,j) * uu(i+1,j) + ss(2,i,j) * uu(i-1,j)
          dd = dd + ss(3,i,j) * uu(i,j+1) + ss(4,i,j) * uu(i,j-1)
          if ( nx > 1 ) then
             if (bc_skewed(mm(i,j),1,+1)) then
                dd = dd + ss(XBC,i,j)*uu(i+2,j)
             else if (bc_skewed(mm(i,j),1,-1)) then
                dd = dd + ss(XBC,i,j)*uu(i-2,j)
             end if
          end if
          if ( ny > 1 ) then
             if (bc_skewed(mm(i,j),2,+1)) then
                dd = dd + ss(YBC,i,j)*uu(i,j+2)
             else if (bc_skewed(mm(i,j),2,-1)) then
                dd = dd + ss(YBC,i,j)*uu(i,j-2)
             end if
          end if
          wrk(i,j) = ff(i,j) - dd
       end do
    end do

    do j = 1, ny
       do i = 1, nx
          uu(i,j) = uu(i,j) + omega*(wrk(i,j)/ss(0,i,j)-uu(i,j))
       end do
    end do

    call destroy(bpt)

  end subroutine jac_smoother_2d

  subroutine jac_dense_smoother_2d(ss, uu, ff, ng)
    use bl_prof_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: ss(0:,:,:)
    real (kind = dp_t), intent(in) :: ff(:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:)
    real (kind = dp_t) :: wrk(size(ff,1),size(ff,2))
    integer :: i, j
    integer :: nx, ny

    real (kind = dp_t), parameter :: omega = 4.0_dp_t/5.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_dense_smoother_2d")

    nx=size(ss,dim=2)
    ny=size(ss,dim=3)

    do j = 1, ny
       do i = 1, nx
          wrk(i,j) = ff(i,j) - (&
               + ss(1,i,j)*uu(i-1,j-1) &
               + ss(2,i,j)*uu(i  ,j-1) &
               + ss(3,i,j)*uu(i+1,j-1) &
               
               + ss(4,i,j)*uu(i-1,j  ) &
               + ss(5,i,j)*uu(i+1,j  ) &
               
               + ss(6,i,j)*uu(i-1,j+1) &
               + ss(7,i,j)*uu(i  ,j+1) &
               + ss(8,i,j)*uu(i+1,j+1) &
               )
       end do
    end do

    do j = 1, ny
       do i = 1, nx
          uu(i,j) = uu(i,j) + omega*(wrk(i,j)/ss(0,i,j)-uu(i,j))
       end do
    end do

    call destroy(bpt)

  end subroutine jac_dense_smoother_2d

  subroutine jac_smoother_3d(ss, uu, ff, mm, ng)
    use bl_prof_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    integer           , intent(in) :: mm(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), allocatable :: wrk(:,:,:)
    real (kind = dp_t) :: dd
    integer :: nx, ny, nz
    integer :: i, j, k
    integer, parameter ::  XBC = 7, YBC = 8, ZBC = 9

    real (kind = dp_t), parameter :: omega = 6.0_dp_t/7.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_smoother_3d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)
    nz = size(ff,dim=3)

    allocate(wrk(nx,ny,nz))

    !$OMP PARALLEL DO PRIVATE(j,i,k,dd) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             dd = ss(1,i,j,k)*uu(i+1,j,k) + ss(2,i,j,k)*uu(i-1,j,k) + &
                  ss(3,i,j,k)*uu(i,j+1,k) + ss(4,i,j,k)*uu(i,j-1,k) + &
                  ss(5,i,j,k)*uu(i,j,k+1) + ss(6,i,j,k)*uu(i,j,k-1)

             if ( nx > 1 ) then
                if (bc_skewed(mm(i,j,k),1,+1)) then
                   dd = dd + ss(XBC,i,j,k)*uu(i+2,j,k)
                else if (bc_skewed(mm(i,j,k),1,-1)) then
                   dd = dd + ss(XBC,i,j,k)*uu(i-2,j,k)
                end if
             end if

             if ( ny > 1 ) then
                if (bc_skewed(mm(i,j,k),2,+1)) then
                   dd = dd + ss(YBC,i,j,k)*uu(i,j+2,k)
                else if (bc_skewed(mm(i,j,k),2,-1)) then
                   dd = dd + ss(YBC,i,j,k)*uu(i,j-2,k)
                end if
             end if

             if ( nz > 1 ) then
                if (bc_skewed(mm(i,j,k),3,+1)) then
                   dd = dd + ss(ZBC,i,j,k)*uu(i,j,k+2)
                else if (bc_skewed(mm(i,j,k),3,-1)) then
                   dd = dd + ss(ZBC,i,j,k)*uu(i,j,k-2)
                end if
             end if
             wrk(i,j,k) = ff(i,j,k) - dd
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             uu(i,j,k) = uu(i,j,k) + omega*(wrk(i,j,k)/ss(0,i,j,k) - uu(i,j,k))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine jac_smoother_3d

  subroutine jac_dense_smoother_3d(ss, uu, ff, ng)
    use bl_prof_module
    integer, intent(in) :: ng
    real (kind = dp_t), intent(in) :: ss(0:,:,:,:)
    real (kind = dp_t), intent(in) :: ff(:,:,:)
    real (kind = dp_t), intent(inout) :: uu(1-ng:,1-ng:,1-ng:)
    real (kind = dp_t), allocatable :: wrk(:,:,:)
    integer :: nx, ny, nz
    integer :: i, j, k

    real (kind = dp_t), parameter :: omega = 6.0_dp_t/7.0_dp_t

    type(bl_prof_timer), save :: bpt

    call build(bpt, "jac_dense_smoother_3d")

    nx = size(ff,dim=1)
    ny = size(ff,dim=2)
    nz = size(ff,dim=3)

    allocate(wrk(nx,ny,nz))

    !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             wrk(i,j,k) = ff(i,j,k) - ( &
                  + ss( 0,i,j,k)*uu(i  ,j  ,k)  &
                  + ss( 1,i,j,k)*uu(i-1,j-1,k-1) &
                  + ss( 2,i,j,k)*uu(i  ,j-1,k-1) &
                  + ss( 3,i,j,k)*uu(i+1,j-1,k-1) &
                  + ss( 4,i,j,k)*uu(i-1,j  ,k-1) &
                  + ss( 5,i,j,k)*uu(i  ,j  ,k-1) &
                  + ss( 6,i,j,k)*uu(i+1,j  ,k-1) &
                  + ss( 7,i,j,k)*uu(i-1,j+1,k-1) &
                  + ss( 8,i,j,k)*uu(i  ,j+1,k-1) &
                  + ss( 9,i,j,k)*uu(i+1,j+1,k-1) &
                  
                  + ss(10,i,j,k)*uu(i-1,j-1,k  ) &
                  + ss(11,i,j,k)*uu(i  ,j-1,k  ) &
                  + ss(12,i,j,k)*uu(i+1,j-1,k  ) &
                  + ss(13,i,j,k)*uu(i-1,j  ,k  ) &
                  + ss(14,i,j,k)*uu(i+1,j  ,k  ) &
                  + ss(15,i,j,k)*uu(i-1,j+1,k  ) &
                  + ss(16,i,j,k)*uu(i  ,j+1,k  ) &
                  + ss(17,i,j,k)*uu(i+1,j+1,k  ) &
                  
                  + ss(18,i,j,k)*uu(i-1,j-1,k+1) &
                  + ss(19,i,j,k)*uu(i  ,j-1,k+1) &
                  + ss(20,i,j,k)*uu(i+1,j-1,k+1) &
                  + ss(21,i,j,k)*uu(i-1,j  ,k+1) &
                  + ss(22,i,j,k)*uu(i  ,j  ,k+1) &
                  + ss(23,i,j,k)*uu(i+1,j  ,k+1) &
                  + ss(24,i,j,k)*uu(i-1,j+1,k+1) &
                  + ss(25,i,j,k)*uu(i  ,j+1,k+1) &
                  + ss(26,i,j,k)*uu(i+1,j+1,k+1) &
                  )
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO PRIVATE(j,i,k) IF(nz.ge.4)
    do k = 1, nz
       do j = 1, ny
          do i = 1, nx
             uu(i,j,k) = uu(i,j,k) + omega*(wrk(i,j,k)/ss(14,i,j,k)-uu(i,j,k))
          end do
       end do
    end do
    !$OMP END PARALLEL DO

    call destroy(bpt)

  end subroutine jac_dense_smoother_3d

end module cc_smoothers_module
