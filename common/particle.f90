module particle
  use cudafor
  implicit none

  private

  public :: particle__init, particle__solv

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  real(8), save :: delx, delt, u0, c, d_delx
  real(8), allocatable :: q(:), r(:)
  
  
contains


  subroutine particle__init(ndim_in,np_in,nsp_in,nxgs_in,nxge_in,nygs_in,nyge_in,nys_in,nye_in, &
                            delx_in,delt_in,c_in,q_in,r_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)

    ndim  = ndim_in
    np    = np_in
    nsp   = nsp_in
    nxgs  = nxgs_in
    nxge  = nxge_in
    nygs  = nygs_in
    nyge  = nyge_in
    nys   = nys_in
    nye   = nye_in
    delx  = delx_in
    delt  = delt_in
    c     = c_in
    allocate(q(nsp))
    allocate(r(nsp))
    q     = q_in
    r     = r_in
    d_delx = 1./delx
  
    is_init = .true.  
  
  end subroutine particle__init


  subroutine particle__solv(gp,up,uf,cumcnt,nxs,nxe)
    implicit none
    integer, intent(in)          :: nxs, nxe
    integer, device, intent(in)  :: cumcnt(nxgs:nxge+1,nys:nye,nsp)
    real(8), device, intent(in)  :: up(ndim,np,nys:nye,nsp)
    real(8), device, intent(in)  :: uf(6,nxgs-2:nxge+2,nys-2:nye+2)
    real(8), device, intent(out) :: gp(ndim,np,nys:nye,nsp)

    integer                      :: i, j, ii, isp
    real(8), device              :: tmp(6,nxs-1:nxe+1,nys-1:nye+1)

    type(dim3) :: Th, Bl

    if(.not.is_init)then
       write(6,*)'Initialize first by calling particle__init()'
       stop
    endif
    
    !fields at (i+1/2, j+1/2)
    !$cuf kernel do (2) <<<*, *>>>
    do j=nys-1,nye+1
    do i=nxs-1,nxe+1
       tmp(1,i,j) = 0.5*(+uf(1,i,j)+uf(1,i,j+1))
       tmp(2,i,j) = 0.5*(+uf(2,i,j)+uf(2,i+1,j))
       tmp(3,i,j) = 0.25*(+uf(3,i,j)  +uf(3,i+1,j) &
                          +uf(3,i,j+1)+uf(3,i+1,j+1))
       tmp(4,i,j) = 0.5*(+uf(4,i,j)+uf(4,i+1,j))
       tmp(5,i,j) = 0.5*(+uf(5,i,j)+uf(5,i,j+1))
       tmp(6,i,j) = uf(6,i,j)
    enddo
    enddo

    Th = dim3(64,1,1)
    Bl = Bl = dim3(ceiling(real(nxe-nxs+1)/Th%x),ceiling(real(nye-nys+1)/Th%y),1)
    call particle__solv_ker<<<Bl,Th>>>(gp,up,uf,tmp,cumcnt,ndim,np,nsp,nxgs,nxge,nxs,nxe,nys,nye,q,r,delt,d_delx)

    if(ndim == 6)then
      !$cuf kernel do (3) <<<*,*>>>
      do isp=1,nsp
        do j = nys, nye
          do ii=1,np
            gp(6,ii,j,isp) = up(6,ii,j,isp)
          enddo
        enddo
      enddo
    endif

  end subroutine particle__solv
  
  attributes(global) &
  subroutine particle__solv_ker(gp,up,uf,tmp,cumcnt,ndim,np,nsp,nsge,nxge,nxs,nxe,nys,nye,q,r,delt,d_delx)
    implicit none
    integer, value :: ndim, np, nsp, nxgs, nxge, nxs, nxe, nys, nye
    real(8), value :: q(nsp), r(nsp), delt, d_delx

    integer, device, intent(in) :: cumcnt(nxgs:nxge+1,nys:nye,nsp)
    real(8), device, intent(in) :: up(ndim, np, nys:nye, nsp)
    real(8), device, intent(in) :: uf(6, nxgs-2:nxge+2, nys-2:nye+2)
    real(8), device, intent(out) :: gp(ndim, np, nys:nye, nsp)

    real(8), device, intent(in) :: tmp(6, nxs-1:nxe+1, nys-1:nye+1)

    integer :: i, j, ii, isp
    real(8) :: dh, sh(-1:1, 2)
    real(8) :: fac1, fac2, txxx, fac1r, fac2r, gam, igam
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6

    i = (blockIdx%x-1)*blockDim%x+threadIdx%x+nxs-1
    j = (blockIdx%y-1)*blockDim%y+threadIdx%y+nys-1

    if (nxs <= i .and. i <= nxe .and. nys <= j .and. j <= nye) then
      do isp=1,nsp

        fac1 = q(isp)/r(isp)*0.5*delt
        txxx = fac1*fac1
        fac2 = q(isp)*delt/r(isp)

        do ii=cumcnt(i,j,isp)+1,cumcnt(i+1,j,isp)

           !second order shape function
           dh = up(1,ii,j,isp)*d_delx-0.5-i
           sh(-1,1) = 0.5*(0.5-dh)*(0.5-dh)
           sh( 0,1) = 0.75-dh*dh
           sh(+1,1) = 0.5*(0.5+dh)*(0.5+dh)

           dh = up(2,ii,j,isp)*d_delx-0.5-j
           sh(-1,2) = 0.5*(0.5-dh)*(0.5-dh)
           sh( 0,2) = 0.75-dh*dh
           sh(+1,2) = 0.5*(0.5+dh)*(0.5+dh)

           bpx = +(+tmp(1,i-1,j-1)*sh(-1,1)+tmp(1,i,j-1)*sh(0,1)+tmp(1,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                 +(+tmp(1,i-1,j  )*sh(-1,1)+tmp(1,i,j  )*sh(0,1)+tmp(1,i+1,j  )*sh(+1,1))*sh( 0,2) &
                 +(+tmp(1,i-1,j+1)*sh(-1,1)+tmp(1,i,j+1)*sh(0,1)+tmp(1,i+1,j+1)*sh(+1,1))*sh(+1,2)

           bpy = +(+tmp(2,i-1,j-1)*sh(-1,1)+tmp(2,i,j-1)*sh(0,1)+tmp(2,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                 +(+tmp(2,i-1,j  )*sh(-1,1)+tmp(2,i,j  )*sh(0,1)+tmp(2,i+1,j  )*sh(+1,1))*sh( 0,2) &
                 +(+tmp(2,i-1,j+1)*sh(-1,1)+tmp(2,i,j+1)*sh(0,1)+tmp(2,i+1,j+1)*sh(+1,1))*sh(+1,2)

           bpz = +(+tmp(3,i-1,j-1)*sh(-1,1)+tmp(3,i,j-1)*sh(0,1)+tmp(3,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                 +(+tmp(3,i-1,j  )*sh(-1,1)+tmp(3,i,j  )*sh(0,1)+tmp(3,i+1,j  )*sh(+1,1))*sh( 0,2) &
                 +(+tmp(3,i-1,j+1)*sh(-1,1)+tmp(3,i,j+1)*sh(0,1)+tmp(3,i+1,j+1)*sh(+1,1))*sh(+1,2)

           epx = +(+tmp(4,i-1,j-1)*sh(-1,1)+tmp(4,i,j-1)*sh(0,1)+tmp(4,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                 +(+tmp(4,i-1,j  )*sh(-1,1)+tmp(4,i,j  )*sh(0,1)+tmp(4,i+1,j  )*sh(+1,1))*sh( 0,2) &
                 +(+tmp(4,i-1,j+1)*sh(-1,1)+tmp(4,i,j+1)*sh(0,1)+tmp(4,i+1,j+1)*sh(+1,1))*sh(+1,2)

           epy = +(+tmp(5,i-1,j-1)*sh(-1,1)+tmp(5,i,j-1)*sh(0,1)+tmp(5,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                 +(+tmp(5,i-1,j  )*sh(-1,1)+tmp(5,i,j  )*sh(0,1)+tmp(5,i+1,j  )*sh(+1,1))*sh( 0,2) &
                 +(+tmp(5,i-1,j+1)*sh(-1,1)+tmp(5,i,j+1)*sh(0,1)+tmp(5,i+1,j+1)*sh(+1,1))*sh(+1,2)

           epz = +(+tmp(6,i-1,j-1)*sh(-1,1)+tmp(6,i,j-1)*sh(0,1)+tmp(6,i+1,j-1)*sh(+1,1))*sh(-1,2) &
                 +(+tmp(6,i-1,j  )*sh(-1,1)+tmp(6,i,j  )*sh(0,1)+tmp(6,i+1,j  )*sh(+1,1))*sh( 0,2) &
                 +(+tmp(6,i-1,j+1)*sh(-1,1)+tmp(6,i,j+1)*sh(0,1)+tmp(6,i+1,j+1)*sh(+1,1))*sh(+1,2)

           !accel.
           uvm1 = up(3,ii,j,isp)+fac1*epx
           uvm2 = up(4,ii,j,isp)+fac1*epy
           uvm3 = up(5,ii,j,isp)+fac1*epz

           !rotate
           gam = dsqrt(c*c+uvm1*uvm1+uvm2*uvm2+uvm3*uvm3)
           igam = 1./gam
           fac1r = fac1*igam
           fac2r = fac2/(gam+txxx*(bpx*bpx+bpy*bpy+bpz*bpz)*igam)

           uvm4 = uvm1+fac1r*(+uvm2*bpz-uvm3*bpy)
           uvm5 = uvm2+fac1r*(+uvm3*bpx-uvm1*bpz)
           uvm6 = uvm3+fac1r*(+uvm1*bpy-uvm2*bpx)

           uvm1 = uvm1+fac2r*(+uvm5*bpz-uvm6*bpy)
           uvm2 = uvm2+fac2r*(+uvm6*bpx-uvm4*bpz)
           uvm3 = uvm3+fac2r*(+uvm4*bpy-uvm5*bpx)

           !accel.
           gp(3,ii,j,isp) = uvm1+fac1*epx
           gp(4,ii,j,isp) = uvm2+fac1*epy
           gp(5,ii,j,isp) = uvm3+fac1*epz

           !move
           gam = 1./dsqrt(1.0+(+gp(3,ii,j,isp)*gp(3,ii,j,isp) &
                                +gp(4,ii,j,isp)*gp(4,ii,j,isp) &
                                +gp(5,ii,j,isp)*gp(5,ii,j,isp))/(c*c))

           gp(1,ii,j,isp) = up(1,ii,j,isp)+gp(3,ii,j,isp)*delt*gam
           gp(2,ii,j,isp) = up(2,ii,j,isp)+gp(4,ii,j,isp)*delt*gam

        enddo
     
     enddo
    end if
      
  end subroutine particle__solv_ker
end module particle
