module sort

  implicit none

  private

  public :: sort__init, sort__bucket

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye


contains


  subroutine sort__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in

    is_init = .true.

  end subroutine sort__init


  subroutine sort__bucket(gp,up,cumcnt,np2,nxs,nxe)
    implicit none
    !BUCKET SORT FOR PARTICLES IN X
    integer, intent(in)          :: nxs, nxe
    integer, device, intent(in)  :: np2(nys:nye,nsp)
    integer, device, intent(out) :: cumcnt(nxgs:nxge+1,nys:nye,nsp)
    real(8), device, intent(in)  :: up(ndim,np,nys:nye,nsp)
    real(8), device, intent(out) :: gp(ndim,np,nys:nye,nsp)
    
    integer                      :: isp
    integer, device              :: cnt(nxs:nxe, nys:nye), sum_cnt(nxs:nxe+1, nys:nye)

    type(dim3)                   :: Th, Bl

    if(.not.is_init)then
       write(6,*)'Initialize first by calling sort__init()'
       stop
    endif

    do isp=1,nsp
      cnt = 0

      Th = dim3(64,1,1)
      Bl = dim3(ceiling(real(nye-nys+1)/Th%x),ceiling(real(np)/Th%y),1)
      call sort__bucket_ker1<<<Bl,Th>>>(up,np2,cnt,ndim,np,nsp,nxs,nxe,nys,nye,isp)

      Th = dim3(64,1,1)
      Bl = dim3(ceiling(real(nye-nys+1)/Th%x),1,1)
      call sort__bucket_ker2<<<Bl,Th>>>(cumcnt,cnt,sum_cnt,nsp,nxgs,nxge,nxs,nxe,nys,nye,isp)
      
      Th = dim3(64,1,1)
      Bl = dim3(ceiling(real(nye-nys+1)/Th%x),ceiling(real(np)/Th%y),1)
      call sort__bucket_ker3<<<Bl,Th>>>(gp,up,np2,sum_cnt,ndim,np,nsp,nxs,nxe,nys,nye,isp)
    enddo

  end subroutine sort__bucket

  attributes(global) &
  subroutine sort__bucket_ker1(up,np2,cnt,ndim,np,nsp,nxs,nxe,nys,nye,isp)
    implicit none
    integer, value                 :: ndim, np, nsp, nxs, nxe, nys, nye, isp
    real(8), device, intent(in)    :: up(ndim,np,nys:nye,nsp)
    integer, device, intent(in)    :: np2(nys:nye,nsp)
    integer, device, intent(inout) :: cnt(nxs:nxe,nys:nye)
    integer                        :: i, j, ii, tmp

    j = (blockIdx%x-1)*blockDim%x+threadIdx%x+nys-1
    ii = (blockIdx%y-1)*blockDim%y+threadIdx%y

    if (nys <= j .and. j <= nye .and. 1 <= ii .and. ii <= np2(j, isp)) then
      i = int(up(1, ii, j, isp))
      tmp = atomicadd(cnt(i, j), 1)
    endif

  end subroutine sort__bucket_ker1

  attributes(global) &
  subroutine sort__bucket_ker2(cumcnt,cnt,sum_cnt,nsp,nxgs,nxge,nxs,nxe,nys,nye,isp)
    implicit none
    integer, value               :: nsp, nxgs, nxge, nxs, nxe, nys, nye, isp
    integer, device, intent(out) :: cumcnt(nxgs:nxge+1,nys:nye,nsp)
    integer, device, intent(in)  :: cnt(nxs:nxe,nys:nye)
    integer, device, intent(out) :: sum_cnt(nxs:nxe+1,nys:nye)
    integer                      :: i, j, ii

    j = (blockIdx%x-1)*blockDim%x+threadIdx%x+nys-1

    if (nys <= j .and. j <= nye) then
      sum_cnt(nxs,j) = 0
      cumcnt(nxs,j,isp) = 0
      do i = nxs+1, nxe+1
        sum_cnt(i,j) = sum_cnt(i-1,j)+cnt(i-1,j)
        cumcnt(i,j,isp) = sum_cnt(i,j)
      enddo
    endif

  end subroutine sort__bucket_ker2

  attributes(global) &
  subroutine sort__bucket_ker3(gp,up,np2,sum_cnt,ndim,np,nsp,nxs,nxe,nys,nye,isp)
    implicit none
    integer, value                 :: ndim, np, nsp, nxs, nxe, nys, nye, isp
    real(8), device, intent(in)    :: up(ndim, np, nys:nye, nsp)
    real(8), device, intent(out)   :: gp(ndim, np, nys:nye, nsp)
    integer, device, intent(in)    :: np2(nys:nye, nsp)
    integer, device, intent(inout) :: sum_cnt(nxs:nxe+1, nys:nye)
    
    integer :: i, j, ii, tmp

    if (nys <= j .and. j <= nye .and. 1 <= ii .and. ii <= np2(j, isp)) then
      i = int(up(1,ii,j,isp))
      tmp = atomicadd(sum_cnt(i,j),1)
      gp(1:ndim,tmp+1,j,isp) = up(1:ndim,ii,j,isp)
    endif

  end subroutine sort__bucket_ker3
end module sort
