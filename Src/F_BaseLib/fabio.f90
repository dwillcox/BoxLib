!! The FABIO_MODULE manages doing input and output of fabs
!! and multifabs in a 'legacy' fashion, that is so that they
!! can be procesed by amrvis?d and the amrderive code.
module fabio_module

  use bl_error_module
  use bl_string_module
  use bl_types
  use fab_module
  use multifab_module
  use parallel

  implicit none

  integer, parameter :: FABIO_DOUBLE = 1
  integer, parameter :: FABIO_SINGLE = 2
  

  interface
     subroutine fabio_close(fd)
       integer, intent(out) :: fd
     end subroutine fabio_close

     subroutine fabio_read_d(fd, offset, d, count)
       use bl_types
       integer, intent(in) :: offset, fd, count
       real(kind=dp_t), intent(out) :: d(count)
     end subroutine fabio_read_d
     subroutine fabio_read_s(fd, offset, s, count)
       use bl_types
       integer, intent(in) :: offset, fd, count
       real(kind=sp_t), intent(out) :: s(count)
     end subroutine fabio_read_s

     subroutine fabio_read_comp_d(fd, offset, skip, d, count)
       use bl_types
       integer, intent(in) :: offset, fd, count, skip
       real(kind=dp_t), intent(out) :: d(count)
     end subroutine fabio_read_comp_d
     subroutine fabio_read_comp_s(fd, offset, skip, s, count)
       use bl_types
       integer, intent(in) :: offset, fd, count, skip
       real(kind=sp_t), intent(out) :: s(count)
     end subroutine fabio_read_comp_s

     subroutine fabio_write_raw_d(fd, offset, d, count, dm, lo, hi, nd, nc, prec)
       use bl_types
       integer, intent(in) :: fd, count, dm, lo(dm), hi(dm), nd(dm), nc, prec
       real(kind=dp_t), intent(in) :: d(count)
       integer, intent(out) :: offset
     end subroutine fabio_write_raw_d
     subroutine fabio_write_raw_s(fd, offset, s, count, dm, lo, hi, nd, nc, prec)
       use bl_types
       integer, intent(in) :: fd, count, dm, lo(dm), hi(dm), nd(dm), nc, prec
       real(kind=sp_t), intent(in) :: s(count)
       integer, intent(out) :: offset
     end subroutine fabio_write_raw_s

  end interface

  private :: fabio_write_raw_d
  private :: fabio_write_raw_s

  interface fabio_write
     module procedure fabio_fab_write_d
     module procedure fabio_multifab_write_d
  end interface

  interface fabio_ml_write
     module procedure fabio_ml_multifab_write_d
  end interface

  integer, parameter :: FABIO_RDONLY = 0
  integer, parameter :: FABIO_WRONLY = 1
  integer, parameter :: FABIO_RDWR   = 2

  integer, parameter :: FABIO_MAX_VAR_NAME = 20
  integer, parameter :: FABIO_MAX_PATH_NAME = 128

contains

  subroutine fabio_mkdir(dirname, stat)
    character(len=*), intent(in) :: dirname
    integer, intent(out), optional :: stat
    interface
       subroutine fabio_mkdir_str(ifilename, stat)
         integer, intent(in) :: ifilename(*)
         integer, intent(out) :: stat
       end subroutine fabio_mkdir_str
    end interface
    integer :: istr(128)
    integer :: lstat

    ! octal conversion 0755
    lstat = 0; if ( present(stat) ) lstat = 1
    call str2int(istr, 128, dirname)
    call fabio_mkdir_str(istr, lstat)
    if ( present(stat) ) stat = lstat

  end subroutine fabio_mkdir

  subroutine fabio_open(fd, filename, mode)
    character(len=*), intent(in):: filename
    integer, intent(out) :: fd
    integer, intent(in), optional :: mode
    interface
       subroutine fabio_open_str(fd, ifilename, mode)
         integer, intent(out) :: fd
         integer, intent(in) :: ifilename(*)
         integer, intent(in) :: mode
       end subroutine fabio_open_str
    end interface
    integer :: istr(128)
    integer :: lmode

    lmode = FABIO_RDONLY
    if ( present(mode) ) then
       if ( mode /= 0 ) lmode = mode
    end if

    call str2int(istr, 128, filename)
    call fabio_open_str(fd, istr, lmode)

  end subroutine fabio_open

  subroutine fabio_fab_write_d(fd, offset, fb, nodal, all, prec)
    integer, intent(in) :: fd
    integer, intent(out) :: offset
    type(fab), intent(in) :: fb
    logical, intent(in), optional :: nodal(:)
    logical, intent(in), optional :: all
    integer, intent(in), optional :: prec
    type(box) :: bx
    logical :: lall
    real(kind=dp_t), pointer :: fbp(:,:,:,:)
    integer :: count, lo(fb%dim), hi(fb%dim), nd(fb%dim), nc, lprec
    lprec = FABIO_DOUBLE; if ( present(prec) ) lprec = prec
    if ( lprec /= FABIO_DOUBLE .and. lprec /= FABIO_SINGLE ) then
       call bl_error("FABIO_WRITE: prec is wrong ", lprec)
    end if
    lall = .false.; if ( present(all) ) lall = all
    bx = get_ibox(fb)
    count = volume(bx)
    fbp => dataptr(fb, bx)
    nc = ncomp(fb)
    lo = lwb(bx)
    hi = upb(bx)
    nd = 0
    if ( present(nodal) ) then
       where ( nodal ) nd = 1
    end if
    call fabio_write_raw_d(fd, offset, fbp, count, fb%dim, lo, hi, nd, nc, lprec)
  end subroutine fabio_fab_write_d

  subroutine fabio_multifab_write_d(mf, dirname, header, prec)
    use bl_IO_module
    type(multifab), intent(in) :: mf
    character(len=*), intent(in) :: dirname, header
    integer, intent(in), optional :: prec
    character(len=128) :: fname
    integer :: un
    integer :: nc, nb, i, fd, j
    integer, allocatable :: offset(:), loffset(:)
    type(box) :: bx
    real(kind=dp_t), allocatable :: mx(:,:), mn(:,:)
    real(kind=dp_t), allocatable :: mxl(:),  mnl(:)
    integer, parameter :: MSG_TAG = 1010

    nc = multifab_ncomp(mf)
    nb = nboxes(mf)
    allocate(offset(nb),loffset(nb))
    allocate(mnl(nc), mxl(nc))
    if ( parallel_IOProcessor() ) then
       allocate(mx(nc,nb), mn(nc,nb))
       un = unit_new()
       call fabio_mkdir(dirname)
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header) // "_H", &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       write(unit=un, fmt='(i0/i0/i0/i0)') 1, 0, nc, 0
       write(unit=un, fmt='("(",i0," 0")')    nb
       do i = 1, nb
          bx = get_box(mf, i)
          call box_print(bx, unit = un, legacy = .True., nodal = mf%nodal)
       end do
       write(unit=un, fmt='(")")')
    end if

    offset = -Huge(offset)
    ! Each processor writes his own FABS
    write(unit=fname, fmt='(a,"_D_",i5.5)') trim(header), parallel_myproc()
    call fabio_open(fd, trim(dirname) // "/" // trim(fname), FABIO_WRONLY)
    do i = 1, nb; if ( remote(mf, i) ) cycle
       call fabio_write(fd, offset(i), mf%fbs(i), nodal = mf%nodal, prec = prec)
    end do
    call fabio_close(fd)
    call parallel_reduce(loffset, offset, MPI_MAX, parallel_IOProcessorNode())

    do i = 1, nb
       if ( local(mf, i) ) then
          do j = 1, nc
             mnl(j) = min_val(mf%fbs(i),j)
             mxl(j) = max_val(mf%fbs(i),j)
          end do
       end if
       if ( parallel_IOProcessor() ) then
          if ( remote(mf, i) ) then
             call parallel_recv(mn(:,i), get_proc(mf%la,i), MSG_TAG)
             call parallel_recv(mx(:,i), get_proc(mf%la,i), MSG_TAG + 1)
          else
             mx(:,i) = mxl
             mn(:,i) = mnl
          end if
       else
          if ( local(mf, i) ) then
            call parallel_send(mnl, parallel_IOProcessorNode(), MSG_TAG)
            call parallel_send(mxl, parallel_IOProcessorNode(), MSG_TAG + 1)
          end if
       end if
    end do

    if ( parallel_IOProcessor() ) then
       if ( any(loffset < 0) ) then
          call bl_error("FABIO_MULTIFAB_WRITE: some loffsets < 0")
       end if
       write(unit=un, fmt='(i0)') nb
       do i = 1, nb
          write(unit=fname, fmt='(a,"_D_",i5.5)') trim(header), get_proc(mf%la, i)
          write(unit=un, fmt='("FabOnDisk: ", a, " ", i0)') trim(fname), loffset(i)
       end do
       write(unit=un, fmt='()')
       write(unit=un, fmt='(i0,",",i0)') nb, nc
       do i = 1, nb
          do j = 1, nc
             write(unit=un, fmt='(es30.20e3,1x,",")') mn(j,i)
          end do
       end do
       write(unit=un, fmt='()')
       write(unit=un, fmt='(i0,",",i0)') nb, nc
       do i = 1, nb
          do j = 1, nc
             write(unit=un, fmt='(es30.20e3,1x,",")') mx(j,i)
          end do
       end do
       close(unit=un)
    end if
  end subroutine fabio_multifab_write_d

  subroutine fabio_ml_multifab_write_d(mfs, rrs, dirname, names, bounding_box, time, dx)
    use bl_IO_module
    type(multifab), intent(in) :: mfs(:)
    integer, intent(in) :: rrs(:)
    character(len=*), intent(in) :: dirname
    type(box), intent(in), optional :: bounding_box
    real(kind=dp_t), intent(in), optional :: time
    real(kind=dp_t), intent(in), optional :: dx(:)
    character(len=FABIO_MAX_VAR_NAME), intent(in), optional :: names(:)
    integer :: i, j, k
    character(len=128) :: header, sd_name
    integer :: nc, un, nl, dm
    real(kind=dp_t), allocatable :: plo(:), phi(:), ldx(:)
    integer, allocatable ::  lo(:),  hi(:)
    integer :: idummy, rdummy
    type(box) :: lbbox
    real(kind=dp_t) :: ltime
    
    if ( size(mfs) < 1 ) then
       call bl_error("FABIO_ML_MULTIFAB_WRITE_D: write a zero length mlmf")
    end if
    if ( size(mfs) /= size(rrs) + 1 ) then
       call bl_error("FABIO_ML_MULTIFAB_WRITE_D: size of mfs /= size (rrs,dim=1)+1")
    end if

    nl = size(mfs)
    nc = ncomp(mfs(1))
    dm = mfs(1)%dim
    allocate(plo(dm),phi(dm),ldx(dm),lo(dm),hi(dm))
    if ( present(bounding_box) ) then
       lbbox = bounding_box
    else
       lbbox = bbox(get_boxarray(mfs(1)))
    end if
    ltime = 0.0_dp_t; if ( present(time) ) ltime = time

    idummy = 0
    rdummy = 0.0_dp_t
    lo = lwb(lbbox); hi = upb(lbbox)
    ldx = 0
    if ( present(dx) ) then
       ldx = dx
    else
       ldx  = 1.0_dp_t/(maxval(hi-lo+1))
    end if
    plo = lwb(lbbox)*ldx
    phi = (upb(lbbox)+1)*ldx

    if ( parallel_IOProcessor() ) then
       call fabio_mkdir(dirname)
    end if

    do i = 1, size(mfs)
       write(unit=sd_name, fmt='(a,"/Level_",i2.2)') trim(dirname), i-1
       call fabio_multifab_write_d(mfs(i), sd_name, "Cell")
    end do

    if ( parallel_IOProcessor() ) then
       header = "Header"
       un = unit_new()
       open(unit=un, &
            file = trim(dirname) // "/" // trim(header), &
            form = "formatted", access = "sequential", &
            status = "replace", action = "write")
       write(unit=un, fmt='("NavierStokes-V1.1")')
       write(unit=un, fmt='(i0)') nc
       if ( present(names) ) then
          do i = 1, nc
             write(unit=un, fmt='(A)') names(i)
          end do
       else
          do i = 1, nc
             write(unit=un, fmt='("Var-",i3.3)') i
          end do
       end if
       write(unit=un, fmt='(i1)') dm
       write(unit=un, fmt='(es30.20e3)') ltime
       write(unit=un, fmt='(i0)') nl - 1
       write(unit=un, fmt='(3es30.20e3)') plo
       write(unit=un, fmt='(3es30.20e3)') phi
       do i = 1, nl - 1
          write(unit=un, fmt='(i0,1x)', advance='no') rrs(i)
       end do
       write(unit=un, fmt='()')
       do i = 1, nl
          call box_print(lbbox, unit=un, legacy = .True., advance = 'no', nodal = mfs(i)%nodal)
          write(unit=un, fmt='()')
          if ( i < nl ) lbbox = refine(lbbox, rrs(i))
       end do
       do i = 1, nl
          write(unit=un, fmt='(i0,1x)', advance = 'no') idummy
       end do
       write(unit=un, fmt='()')
       do i = 1, nl
          write(unit=un, fmt='(3es30.20e3)') ldx
          if ( i < nl ) ldx = ldx/rrs(i)
       end do
       write(unit=un, fmt='(i0)') idummy
       write(unit=un, fmt='(i0)') idummy
       ! SOME STUFF
       do i = 1, nl
          write(unit=un, fmt=*) i-1, nboxes(mfs(i)), rdummy, idummy
          do j = 1, nboxes(mfs(i))
             plo =  lwb(get_box(mfs(i),j))    
             phi = (upb(get_box(mfs(i),j))+1)
             do k = 1, dm
                write(unit=un, fmt='(2es30.20e3)') plo(k)*ldx(k), phi(k)*ldx(k)
             end do
          end do
          write(unit=un, fmt='("Level_",i2.2,"/Cell")') i-1
       end do
       close(unit=un)
    end if
  end subroutine fabio_ml_multifab_write_d

  subroutine fabio_ml_multifab_read_d(mmf, root, unit, ng)
    use bl_stream_module
    use bl_IO_module
    type(multifab), pointer :: mmf(:)
    character(len=*), intent(in) :: root
    integer, intent(in), optional :: unit
    integer, intent(in), optional :: ng
    integer :: lun
    type(bl_stream) :: strm
    character(len=256) :: str

    call build(strm, unit)
    lun = bl_stream_the_unit(strm)
    open(unit=lun, &
         file = trim(root) // "/" // "Header", &
         status = 'old', action = 'read')
    read(unit=lun,fmt='(a)') str
    if ( str == '&ML_MULTIFAB' ) then
       call bl_error("PLOTFILE_BUILD: not implemented")
    else if ( str == 'NavierStokes-V1.1' .or. str == 'HyperCLaw-V1.1' ) then 
       call build_ns_plotfile
    else
       call bl_error('FABIO_ML_MULTIFAB_WRITE_D: Header has improper magic string', str)
    end if
    call destroy(strm)

  contains

    ! NavierStokes-V1.1 Plotfile Formats
    ! Record
    !     : c : NavierStokes-V1.1/HyperClaw-V1.1
    !     : c : Numbers of fields = n
    !    n: i : Field Names
    !     : i : Dimension = dm
    !     : r : Time
    !     : i : Number of Levels - 1 : nl
    !     : r : Physical domain lo end [1:dm]
    !     : r : Physical domain hi end [1:dm]
    !     : i : Refinement Ratios [1:nl-1]
    !     : b : Prob domains per level [1:nl]
    !     : i : unused [1:nl]
    !   nl: r : grid spacing, per level, [1:dm]
    !     : i : unused  :
    !     : i : unused
    !     For each level
    !     [
    !       : iiri : dummy, nboxes, dummy, dummy
    !       For each box, j
    !       [
    !         : r :  plo[1:dm,j], phi[1:dm, j]
    !       ]
    !       : c : level directory
    !     ]
    !     Close Header File
    !     For each level
    !     [
    !       Open Header of sub-directory
    !       : iiii: dummy, dummy, ncomponents, dummy
    !       : i ; '(', nboxes dummy
    !       For each box, j
    !       [
    !         : b : bx[j]
    !       ]
    !       :  : ')'
    !       For each box, j
    !       [
    !         : ci : 'FabOnDisk: ' Filename[j], Offset[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : min[j]
    !       ]
    !       : i : nboxes, ncomponents
    !       For each box, j
    !       [
    !         : r : man[j]
    !       ]
    !       Close subgrid file
    !     ]

    subroutine build_ns_plotfile()
      integer :: i, n, k
      integer :: j, nc
      integer n1
      character(len=FABIO_MAX_PATH_NAME) :: str, str1, cdummy, filename
      integer :: offset
      integer :: idummy, sz, fd
      real(kind=dp_t) :: rdummy, tm
      integer :: nvars, dm, flevel
      integer, allocatable :: refrat(:,:), nboxes(:)
      real(kind=dp_t), allocatable :: dxlev(:,:)
      type(box), allocatable :: pdbx(:), bxs(:)
      type(box) :: bx
      type(boxarray) :: ba
      type(layout) :: la
      character(len=256), allocatable :: fileprefix(:)
      character(len=256), allocatable :: header(:)
      real(kind=dp_t), pointer :: pp(:,:,:,:)

      read(unit=lun,fmt=*) nvars
      do i = 1, nvars
         read(unit=lun,fmt='(a)') cdummy
      end do
      read(unit=lun, fmt=*) dm
      read(unit=lun, fmt=*) tm
      read(unit=lun, fmt=*) flevel
      flevel = flevel + 1

      allocate(mmf(flevel))

      allocate(pdbx(flevel))
      allocate(nboxes(flevel))
      allocate(fileprefix(flevel))
      allocate(header(flevel))

      read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
      !! Not make this really work correctly, I need to see if these are
      !! IntVects here.  I have no examples of this.
      allocate(refrat(flevel-1,1:dm))
      read(unit=lun, fmt=*) refrat(:,1)
      refrat(:,2:dm) = spread(refrat(:,1), dim=2, ncopies=dm-1)

      do i = 1, flevel
         call box_read(pdbx(i), unit = lun)
      end do
      read(unit=lun, fmt=*) (idummy, i=1, flevel)
      allocate(dxlev(flevel,1:dm))
      do i = 1, flevel
         read(unit=lun, fmt=*) dxlev(i,:)
      end do

      read(unit=lun, fmt=*) idummy, idummy
      do i = 1, flevel
         read(unit=lun, fmt=*) idummy, nboxes(i), rdummy, idummy
         do j = 1, nboxes(i)
            read(unit=lun, fmt=*) (rdummy, k=1, 2*dm)
         end do
         read(unit=lun, fmt='(a)') str
         str1 = str(:index(str, "/")-1)
         fileprefix(i) = str1
         str1 = trim(str(index(str, "/")+1:)) // "_H"
         header(i) = trim(str1)
      end do
      close(unit=lun)
      do i = 1, flevel
         open(unit=lun, &
              action = 'read', &
              status = 'old', file = trim(trim(root) // "/" //  &
              trim(fileprefix(i)) // "/" // &
              trim(header(i))) )
         read(unit=lun, fmt=*) idummy, idummy, nc, idummy
         allocate(bxs(nboxes(i)))
         if ( nc /= nvars ) &
              call bl_error("BUILD_PLOTFILE: unexpected nc", nc)
         call bl_stream_expect(strm, '(')
         n = bl_stream_scan_int(strm)
         if ( n /= nboxes(i) ) &
              call bl_error("BUILD_PLOTFILE: unexpected n", n)
         idummy = bl_stream_scan_int(strm)
         do j = 1, nboxes(i)
            call box_read(bxs(j), unit = lun)
         end do
         call bl_stream_expect(strm, ')')
         call build(ba, bxs)
         call build(la, ba)
         call multifab_build(mmf(i),la,nc = nvars, ng = ng)
         read(unit=lun, fmt=*) idummy
         do j = 1, nboxes(i)
            read(unit=lun, fmt=*) cdummy, &
                 filename, offset
            call fabio_open(fd,                         &
               trim(root) // "/" //                &
               trim(fileprefix(i)) // "/" // &
               trim(filename))
            bx = get_box(mmf(i), j)
            pp => dataptr(mmf(i), j, bx)
            sz = volume(bxs(j))
            call fabio_read_d(fd, offset, pp(:,:,:,:), sz*nvars)
            call fabio_close(fd)
         end do
         deallocate(bxs)
         call destroy(ba)
         close(unit=lun)
      end do
    end subroutine build_ns_plotfile
  end subroutine fabio_ml_multifab_read_d

end module fabio_module
