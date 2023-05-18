module boundary_periodic
  use mpi
  use cudafor
  implicit none

  private

  public :: boundary_periodic__init
  public :: boundary_periodic__dfield
  public :: boundary_periodic__particle_x
  public :: boundary_periodic__particle_y
  public :: boundary_periodic__curre
  public :: boundary_periodic__phi
  public :: boundary_periodic__mom

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  integer, save :: nup, ndown, mnpi, mnpr, ncomw
  integer       :: nerr
  integer, allocatable :: nstat(:)
  real(8), save :: delx, delt, c


contains


  subroutine boundary_periodic__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
                            nup_in,ndown_in,mnpi_in,mnpr_in,ncomw_in,nerr_in,nstat_in,          &
                            delx_in,delt_in,c_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nup_in, ndown_in, mnpi_in, mnpr_in, ncomw_in, nerr_in, nstat_in(:)
    real(8), intent(in) :: delx_in, delt_in, c_in

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in
    nup   = nup_in
    ndown = ndown_in
    mnpi  = mnpi_in
    mnpr  = mnpr_in
    ncomw = ncomw_in
    nerr  = nerr_in
    allocate(nstat(size(nstat_in)))
    nstat = nstat_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in

    is_init = .true.

  end subroutine boundary_periodic__init


  subroutine boundary_periodic__particle_x(up,np2)
    implicit none
    integer, device, intent(in)    :: np2(nys:nye,nsp)
    real(8), device, intent(inout) :: up(ndim,np,nys:nye,nsp)
    type(dim3)                     :: Th, Bl
    integer                :: j, ii, isp, ipos

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary_periodic__init()'
       stop
    endif

    Th = dim3(64,1,1)
    Bl = dim3(ceiling(real(nye-nys+1)/Th%x), ceiling(real(np)/Th%y), 1)
    do isp=1,nsp
      call boundary_periodic__particle_x_ker<<<Bl,Th>>>(up,np2,ndim,np,nsp,nxgs,nxge,nys,nye,isp,delx)
    enddo

  end subroutine boundary_periodic__particle_x


  attributes(global) &
  subroutine boundary_periodic__particle_x_ker(up,np2,ndim,np,nsp,nxgs,nxge,nys,nye,isp,delx)
    use cudadevice
    implicit none
    integer, value                 :: ndim, np, nsp, nxgs, nxge, nys, nye, isp
    real(8), value                 :: delx
    integer, device, intent(in)    :: np2(nys:nye,nsp)
    real(8), device, intent(inout) :: up(ndim,np,nys:nye,nsp)
    integer                        :: j, ii, ipos

    j = (blockIdx%x-1)*blockDim%x+threadIdx%x+nys-1
    ii = (blockIdx%y-1)*blockDim%y+threadIdx%y

    if(nys <= j .and. j <= nye .and. 1 <= ii .and. ii <= np2(j, isp))then
      ipos = int(__ddiv_rd(up(1, ii, j, isp),delx))

      if(ipos < nxgs)then
        up(1, ii, j, isp) = __fma_rd(dble(nxge-nxgs+1),delx,up(1,ii,j,isp))
      else if(ipos >= nxge+1)then
        up(1, ii, j, isp) = __fma_rd(-dble(nxge-nxgs+1),delx,up(1,ii,j,isp))
      endif
    endif

  end subroutine boundary_periodic__particle_x_ker


  subroutine boundary_periodic__particle_y(up,np2)
    implicit none
    integer, device, intent(inout)     :: np2(nys:nye,nsp)
    real(8), device, intent(inout)     :: up(ndim,np,nys:nye,nsp)
    logical, save                      :: lflag=.true.
    integer                            :: j, isp, cnt_tmp_h, cnt1_h, cnt2_h
    integer, device                    :: cnt(nys-1:nye+1), cnt2(nys:nye), cnt_tmp
    integer, device, save, allocatable :: flag(:,:)
    real(8), device, save, allocatable :: bff_ptcl(:,:)
    type(dim3)                         :: Th, Bl

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary_periodic__init()'
       stop
    endif

    if(lflag)then
       allocate(flag(np,nys:nye))
       allocate(bff_ptcl(ndim*np,nys-1:nye+1))
       lflag=.false.
    endif

    do isp=1,nsp
       cnt = 0
       cnt2 = 0

       Th = dim3(64,1,1)
       Bl = dim3(ceiling(real(nye-nys+1)/Th%x),1,1)
       call boundary_periodic__particle_y_ker1<<<Bl,Th>>>(up,np2,cnt,cnt2,flag,bff_ptcl,ndim,np,nsp,nys,nye,isp,delx)

       !transfer to rank-1
       call MPI_SENDRECV(cnt(nys-1),1,mnpi,ndown,100, &
                         cnt_tmp   ,1,mnpi,nup  ,100, &
                         ncomw,nstat,nerr)

       cnt_tmp_h = cnt_tmp
       cnt1_h = cnt(nys-1)
       cnt2_h = cnt(nye)

       call MPI_SENDRECV(bff_ptcl(1          ,nys-1),ndim*cnt1_h   ,mnpr,ndown,101, &
                         bff_ptcl(ndim*cnt2_h+1,nye),ndim*cnt_tmp_h,mnpr,nup  ,101, &
                         ncomw,nstat,nerr)

       cnt_tmp_h = cnt_tmp_h+cnt(nye)
       cnt(nye) = cnt_tmp_h

       !transfer to rank+1
       call MPI_SENDRECV(cnt(nye+1),1,mnpi,nup  ,200, &
                         cnt_tmp   ,1,mnpi,ndown,200, &
                         ncomw,nstat,nerr)

       cnt_tmp_h = cnt_tmp
       cnt1_h = cnt(nye+1)
       cnt2_h = cnt(nys)

       call MPI_SENDRECV(bff_ptcl(1          ,nye+1),ndim*cnt1_h   ,mnpr,nup  ,201, &
                         bff_ptcl(ndim*cnt2_h+1,nys),ndim*cnt_tmp_h,mnpr,ndown,201, &
                         ncomw,nstat,nerr)

       cnt_tmp_h = cnt_tmp_h+cnt(nys)
       cnt(nys) = cnt_tmp_h

       Th = dim3(64,1,1)
       Bl = dim3(ceiling(real(nye-nys+1)/Th%x),1,1)
       call boundary_periodic__particle_y_ker2<<<Bl,Th>>>(up,np2,cnt,cnt2,flag,bff_ptcl,ndim,np,nsp,nys,nye,isp)

       !$cuf kernel do <<<*,*>>>
       do j=nys,nye
          np2(j,isp) = np2(j,isp)+cnt(j)
          if(np2(j,isp) > np) then
             write(*,*)"memory over (np2 > np)",np,np2(j,isp),j,isp
             stop
          endif
       enddo

    enddo

  end subroutine boundary_periodic__particle_y


  attributes(global) &
  subroutine boundary_periodic__particle_y_ker1(up,np2,cnt,cnt2,flag,bff_ptcl,ndim,np,nsp,nys,nye,isp,delx)
    use cudadevice
    implicit none
    integer, value                 :: ndim, np, nsp, nys, nye, isp
    real(8), value                 :: delx
    integer, device, intent(in)    :: np2(nys:nye, nsp)
    real(8), device, intent(inout) :: up(ndim,np,nys:nye,nsp)
    integer, device, intent(inout) :: cnt(nys-1:nye+1), cnt2(nys:nye)
    integer, device, intent(inout) :: flag(np,nys:nye)
    real(8), device, intent(inout) :: bff_ptcl(ndim*np,nys-1:nye+1)
    integer                        :: tmp, j, ii, jpos, idim

    j = (blockIdx%x-1)*blockDim%x+threadIdx%x+nys-1
    if (nys <= j .and. j <= nye) then
      do ii = 1, np2(j,isp)
        jpos = int(__ddiv_rd(up(2,ii,j,isp),delx))

        if (jpos /= j) then
          if (jpos <= nys-1) then
            up(2,ii,j,isp) = __fma_rd(dble(nye-nys+1),delx,up(2,ii,j,isp))
          else if (jpos >= nye+1) then
            up(2,ii,j,isp) = __fma_rd(-dble(nye-nys+1),delx,up(2,ii,j,isp))
          endif

          tmp = atomicadd(cnt(jpos),1)
          do idim = 1, ndim
            bff_ptcl(idim+ndim*tmp,jpos) = up(idim,ii,j,isp)
          enddo

          tmp = atomicadd(cnt2(j),1)
          flag(tmp+1, j) = ii
        endif
      enddo
    endif

  end subroutine boundary_periodic__particle_y_ker1


  attributes(global) &
  subroutine boundary_periodic__particle_y_ker2(up,np2,cnt,cnt2,flag,bff_ptcl,ndim,np,nsp,nys,nye,isp)
    implicit none
    integer, value                 :: ndim, np, nsp, nys, nye, isp
    integer, device, intent(inout) :: np2(nys:nye,nsp)
    real(8), device, intent(inout) :: up(ndim,np,nys:nye,nsp)
    integer, device, intent(inout) :: cnt(nys-1:nye+1)
    integer, device, intent(in)    :: cnt2(nys:nye)
    integer, device, intent(in)    :: flag(np, nys:nye)
    real(8), device, intent(in)    :: bff_ptcl(ndim*np,nys-1:nye+1)
    integer                        :: cnt_tmp, idim, iii, ii, j

    j = (blockIdx%x-1)*blockDim%x+threadIdx%x+nys-1
    if(nys <= j .and. j <= nye)then
      iii = 0
      cnt_tmp = cnt2(j)
      loop1: do ii = 1,cnt2(j)
        if(cnt(j) == 0)then
          if (np2(j,isp) < flag(ii,j)) exit loop1
          do while (np2(j,isp) == flag(cnt_tmp,j))
            np2(j,isp) = np2(j,isp)-1
            if (np2(j,isp) < flag(ii,j)) exit loop1
            cnt_tmp = cnt_tmp-1
          end do
          do idim = 1,ndim
            up(idim,flag(ii,j),j,isp) = up(idim, np2(j,isp),j,isp)
          end do
          np2(j, isp) = np2(j, isp)-1
        else
          do idim = 1, ndim
            up(idim, flag(ii, j), j, isp) = bff_ptcl(idim+ndim*iii,j)
          end do
          iii = iii+1
          cnt(j) = cnt(j)-1
        end if
      enddo loop1

      if(cnt(j) > 0)then
        do ii = 1,cnt(j)
          do idim = 1,ndim
            up(idim, np2(j,isp)+ii,j,isp) = bff_ptcl(ndim*iii+idim+ndim*(ii-1),j)
          enddo
        enddo
      endif
    endif


  end subroutine boundary_periodic__particle_y_ker2


  subroutine boundary_periodic__dfield(df,nxs,nxe,nys,nye,nxgs,nxge)

    integer, intent(in)            :: nxs, nxe, nys, nye, nxgs, nxge
    real(8), device, intent(inout) :: df(6,nxgs-2:nxge+2,nys-2:nye+2)
    integer                        :: i, j, ii, ieq
    real(8), device                :: bff_snd(12*(nxe-nxs+1)), bff_rcv(12*(nxe-nxs+1))

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary_periodic__init()'
       stop
    endif

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = 12*(i-nxs)
       bff_snd(ii+1)  = df(1,i,nys)
       bff_snd(ii+2)  = df(2,i,nys)
       bff_snd(ii+3)  = df(3,i,nys)
       bff_snd(ii+4)  = df(4,i,nys)
       bff_snd(ii+5)  = df(5,i,nys)
       bff_snd(ii+6)  = df(6,i,nys)
       bff_snd(ii+7)  = df(1,i,nys+1)
       bff_snd(ii+8)  = df(2,i,nys+1)
       bff_snd(ii+9)  = df(3,i,nys+1)
       bff_snd(ii+10) = df(4,i,nys+1)
       bff_snd(ii+11) = df(5,i,nys+1)
       bff_snd(ii+12) = df(6,i,nys+1)
    enddo

    call MPI_SENDRECV(bff_snd(1),12*(nxe-nxs+1),mnpr,ndown,110, &
                      bff_rcv(1),12*(nxe-nxs+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = 12*(i-nxs)
       df(1,i,nye+1) = bff_rcv(ii+1)
       df(2,i,nye+1) = bff_rcv(ii+2)
       df(3,i,nye+1) = bff_rcv(ii+3)
       df(4,i,nye+1) = bff_rcv(ii+4)
       df(5,i,nye+1) = bff_rcv(ii+5)
       df(6,i,nye+1) = bff_rcv(ii+6)
       df(1,i,nye+2) = bff_rcv(ii+7)
       df(2,i,nye+2) = bff_rcv(ii+8)
       df(3,i,nye+2) = bff_rcv(ii+9)
       df(4,i,nye+2) = bff_rcv(ii+10)
       df(5,i,nye+2) = bff_rcv(ii+11)
       df(6,i,nye+2) = bff_rcv(ii+12)
    enddo

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = 12*(i-nxs)
       bff_snd(ii+1)  = df(1,i,nye-1)
       bff_snd(ii+2)  = df(2,i,nye-1)
       bff_snd(ii+3)  = df(3,i,nye-1)
       bff_snd(ii+4)  = df(4,i,nye-1)
       bff_snd(ii+5)  = df(5,i,nye-1)
       bff_snd(ii+6)  = df(6,i,nye-1)
       bff_snd(ii+7)  = df(1,i,nye)
       bff_snd(ii+8)  = df(2,i,nye)
       bff_snd(ii+9)  = df(3,i,nye)
       bff_snd(ii+10) = df(4,i,nye)
       bff_snd(ii+11) = df(5,i,nye)
       bff_snd(ii+12) = df(6,i,nye)
    enddo

    call MPI_SENDRECV(bff_snd(1),12*(nxe-nxs+1),mnpr,nup  ,100, &
                      bff_rcv(1),12*(nxe-nxs+1),mnpr,ndown,100, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = 12*(i-nxs)
       df(1,i,nys-2) = bff_rcv(ii+1)
       df(2,i,nys-2) = bff_rcv(ii+2)
       df(3,i,nys-2) = bff_rcv(ii+3)
       df(4,i,nys-2) = bff_rcv(ii+4)
       df(5,i,nys-2) = bff_rcv(ii+5)
       df(6,i,nys-2) = bff_rcv(ii+6)
       df(1,i,nys-1) = bff_rcv(ii+7)
       df(2,i,nys-1) = bff_rcv(ii+8)
       df(3,i,nys-1) = bff_rcv(ii+9)
       df(4,i,nys-1) = bff_rcv(ii+10)
       df(5,i,nys-1) = bff_rcv(ii+11)
       df(6,i,nys-1) = bff_rcv(ii+12)
    enddo

    !$cuf kernel do (2) <<<*,*>>>
    do j = nys-2,nye+2
      do ieq = 1,6
         df(ieq,nxs-1,j) = df(ieq,nxe,j)
         df(ieq,nxe+1,j) = df(ieq,nxs,j)
      enddo
    enddo

  end subroutine boundary_periodic__dfield


  subroutine boundary_periodic__curre(uj,nxs,nxe,nys,nye,nxgs,nxge)

    integer, intent(in)            :: nxs, nxe, nys, nye, nxgs, nxge
    real(8), device, intent(inout) :: uj(3,nxgs-2:nxge+2,nys-2:nye+2)
    integer                        :: i, ii, j, ieq
    real(8), device                :: bff_rcv(6*(nxe-nxs+4+1)), bff_snd(6*(nxe-nxs+4+1))

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary_periodic__init()'
       stop
    endif

    !send to rank-1
    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nys-2)
       bff_snd(ii+2) = uj(2,i,nys-2)
       bff_snd(ii+3) = uj(3,i,nys-2)
       bff_snd(ii+4) = uj(1,i,nys-1)
       bff_snd(ii+5) = uj(2,i,nys-1)
       bff_snd(ii+6) = uj(3,i,nys-1)
    enddo

    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,ndown,110, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,nup  ,110, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nye-1) = uj(1,i,nye-1)+bff_rcv(ii+1)
       uj(2,i,nye-1) = uj(2,i,nye-1)+bff_rcv(ii+2)
       uj(3,i,nye-1) = uj(3,i,nye-1)+bff_rcv(ii+3)
       uj(1,i,nye  ) = uj(1,i,nye  )+bff_rcv(ii+4)
       uj(2,i,nye  ) = uj(2,i,nye  )+bff_rcv(ii+5)
       uj(3,i,nye  ) = uj(3,i,nye  )+bff_rcv(ii+6)
    enddo

    !send to rank+1
    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nye+1)
       bff_snd(ii+2) = uj(2,i,nye+1)
       bff_snd(ii+3) = uj(3,i,nye+1)
       bff_snd(ii+4) = uj(1,i,nye+2)
       bff_snd(ii+5) = uj(2,i,nye+2)
       bff_snd(ii+6) = uj(3,i,nye+2)
    enddo

    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,nup  ,120, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,ndown,120, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nys  ) = uj(1,i,nys  )+bff_rcv(ii+1)
       uj(2,i,nys  ) = uj(2,i,nys  )+bff_rcv(ii+2)
       uj(3,i,nys  ) = uj(3,i,nys  )+bff_rcv(ii+3)
       uj(1,i,nys+1) = uj(1,i,nys+1)+bff_rcv(ii+4)
       uj(2,i,nys+1) = uj(2,i,nys+1)+bff_rcv(ii+5)
       uj(3,i,nys+1) = uj(3,i,nys+1)+bff_rcv(ii+6)
    enddo

!#####    !Update of nori-shiro   #####

    !send to rank-1
    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nys)
       bff_snd(ii+2) = uj(2,i,nys)
       bff_snd(ii+3) = uj(3,i,nys)
       bff_snd(ii+4) = uj(1,i,nys+1)
       bff_snd(ii+5) = uj(2,i,nys+1)
       bff_snd(ii+6) = uj(3,i,nys+1)
    enddo

    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,ndown,130, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,nup  ,130, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nye+1) = bff_rcv(ii+1)
       uj(2,i,nye+1) = bff_rcv(ii+2)
       uj(3,i,nye+1) = bff_rcv(ii+3)
       uj(1,i,nye+2) = bff_rcv(ii+4)
       uj(2,i,nye+2) = bff_rcv(ii+5)
       uj(3,i,nye+2) = bff_rcv(ii+6)
    enddo

    !send to rank+1
    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       bff_snd(ii+1) = uj(1,i,nye-1)
       bff_snd(ii+2) = uj(2,i,nye-1)
       bff_snd(ii+3) = uj(3,i,nye-1)
       bff_snd(ii+4) = uj(1,i,nye)
       bff_snd(ii+5) = uj(2,i,nye)
       bff_snd(ii+6) = uj(3,i,nye)
    enddo

    call MPI_SENDRECV(bff_snd(1),6*(nxe-nxs+4+1),mnpr,nup  ,140, &
                      bff_rcv(1),6*(nxe-nxs+4+1),mnpr,ndown,140, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs-2,nxe+2
       ii = 6*(i-(nxs-2))
       uj(1,i,nys-2) = bff_rcv(ii+1)
       uj(2,i,nys-2) = bff_rcv(ii+2)
       uj(3,i,nys-2) = bff_rcv(ii+3)
       uj(1,i,nys-1) = bff_rcv(ii+4)
       uj(2,i,nys-1) = bff_rcv(ii+5)
       uj(3,i,nys-1) = bff_rcv(ii+6)
    enddo

    !$cuf kernel do (2) <<<*,*>>>
    do j = nys-2,nye+2
      do ieq = 1,3
         uj(ieq,nxe-1,j) = uj(ieq,nxe-1,j)+uj(ieq,nxs-2,j)
         uj(ieq,nxe  ,j) = uj(ieq,nxe  ,j)+uj(ieq,nxs-1,j)
         uj(ieq,nxs  ,j) = uj(ieq,nxs  ,j)+uj(ieq,nxe+1,j)
         uj(ieq,nxs+1,j) = uj(ieq,nxs+1,j)+uj(ieq,nxe+2,j)
      enddo
   enddo

   !$cuf kernel do (2) <<<*,*>>>
   do j = nys-2,nye+2
      do ieq = 1,3
         uj(ieq,nxs-2,j) = uj(ieq,nxe-1,j)
         uj(ieq,nxs-1,j) = uj(ieq,nxe  ,j)
         uj(ieq,nxe+1,j) = uj(ieq,nxs  ,j)
         uj(ieq,nxe+2,j) = uj(ieq,nxs+1,j)
      enddo
   enddo

  end subroutine boundary_periodic__curre


  subroutine boundary_periodic__phi(phi,nxs,nxe,nys,nye,l)

    integer, intent(in)            :: nxs, nxe, nys, nye, l
    real(8), device, intent(inout) :: phi(nxs-1:nxe+1,nys-1:nye+1)
    integer                        :: i, ii, j
    real(8), device                :: bff_snd(nxe-nxs+1), bff_rcv(nxe-nxs+1)

    if(.not.is_init)then
       write(6,*)'Initialize first by calling boundary_periodic__init()'
       stop
    endif

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = i-nxs
       bff_snd(ii+1)  = phi(i,nys)
    enddo

    call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,ndown,110, &
                      bff_rcv(1),nxe-nxs+1,mnpr,nup  ,110, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = i-nxs
       phi(i,nye+1) = bff_rcv(ii+1)
    enddo

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = i-nxs
       bff_snd(ii+1)  = phi(i,nye)
    enddo

    call MPI_SENDRECV(bff_snd(1),nxe-nxs+1,mnpr,nup  ,100, &
                      bff_rcv(1),nxe-nxs+1,mnpr,ndown,100, &
                      ncomw,nstat,nerr)

    !$cuf kernel do <<<*,*>>>
    do i=nxs,nxe
       ii = i-nxs
       phi(i,nys-1) = bff_rcv(ii+1)
    enddo

    !$cuf kernel do <<<*,*>>>
    do j=nys-1,nye+1
      phi(nxs-1,j) = phi(nxe,j)
      phi(nxe+1,j) = phi(nxs,j)
    enddo

  end subroutine boundary_periodic__phi


  subroutine boundary_periodic__mom(mom)

    real(8), intent(inout) :: mom(7,nxgs-1:nxge+1,nys-1:nye+1,nsp)

    integer, parameter :: nk = 7
    integer :: i, ii, k, isp
    real(8) :: bff_rcv(nk*(nxge-nxgs+3)), bff_snd(nk*(nxge-nxgs+3))

!$OMP PARALLEL WORKSHARE
    mom(1:nk,nxgs,nys-1:nye+1,1:nsp) = mom(1:nk,nxgs  ,nys-1:nye+1,1:nsp) &
                                      +mom(1:nk,nxge+1,nys-1:nye+1,1:nsp)
    mom(1:nk,nxge,nys-1:nye+1,1:nsp) = mom(1:nk,nxge  ,nys-1:nye+1,1:nsp) &
                                      +mom(1:nk,nxgs-1,nys-1:nye+1,1:nsp)
!$OMP END PARALLEL WORKSHARE

    do isp = 1,nsp
      !send to rank-1
!$OMP PARALLEL DO PRIVATE(i,ii)
      do i=nxgs-1,nxge+1
         ii = nk*(i-(nxgs-1))
         do k = 1, nk
            bff_snd(ii+k) = mom(k,i,nys-1,isp)
         enddo
      enddo
!$OMP END PARALLEL DO

      call MPI_SENDRECV(bff_snd(1),nk*(nxge-nxgs+3),mnpr,ndown,200, &
                        bff_rcv(1),nk*(nxge-nxgs+3),mnpr,nup  ,200, &
                        ncomw,nstat,nerr)

!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii)
      do i=nxgs-1,nxge+1
         ii = nk*(i-(nxgs-1))
         do k = 1, nk
            mom(k,i,nye,isp) = mom(k,i,nye,isp)+bff_rcv(ii+k)
         enddo
      enddo
!$OMP END DO NOWAIT

      !send to rank+1
!$OMP DO PRIVATE(i,ii)
      do i=nxgs-1,nxge+1
         ii = nk*(i-(nxgs-1))
         do k = 1, nk
            bff_snd(ii+k) = mom(k,i,nye+1,isp)
         enddo
      enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      call MPI_SENDRECV(bff_snd(1),nk*(nxge-nxgs+3),mnpr,nup  ,201, &
                        bff_rcv(1),nk*(nxge-nxgs+3),mnpr,ndown,201, &
                        ncomw,nstat,nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
      do i=nxgs-1,nxge+1
         ii = nk*(i-(nxgs-1))
         do k = 1, nk
            mom(k,i,nys,isp) = mom(k,i,nys,isp)+bff_rcv(ii+k)
         end do
      enddo
!$OMP END PARALLEL DO
    enddo

  end subroutine boundary_periodic__mom


end module boundary_periodic
