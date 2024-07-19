module particle

  implicit none

  private

  public :: particle__init, particle__solv

  logical, save :: is_init = .false.
  integer, save :: ndim, np, nsp, nxgs, nxge, nygs, nyge, nys, nye
  real(8), save :: delx, delt, u0, c, d_delx
  real(8), allocatable :: q(:), r(:)

contains

  subroutine particle__init(ndim_in, np_in, nsp_in, nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in, &
                            delx_in, delt_in, c_in, q_in, r_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    real(8), intent(in) :: delx_in, delt_in, c_in, q_in(nsp_in), r_in(nsp_in)

    ndim = ndim_in
    np = np_in
    nsp = nsp_in
    nxgs = nxgs_in
    nxge = nxge_in
    nygs = nygs_in
    nyge = nyge_in
    nys = nys_in
    nye = nye_in
    delx = delx_in
    delt = delt_in
    c = c_in
    allocate (q(nsp))
    allocate (r(nsp))
    q = q_in
    r = r_in
    d_delx = 1./delx

    is_init = .true.

  end subroutine particle__init

  subroutine particle__solv(gp, up, uf, cumcnt, nxs, nxe)

    integer, intent(in) :: nxs, nxe
    integer, intent(in) :: cumcnt(nxgs:nxge + 1, nys:nye, nsp)
    real(8), intent(in) :: up(ndim, np, nys:nye, nsp)
    real(8), intent(in) :: uf(6, nxgs - 2:nxge + 2, nys - 2:nye + 2)
    real(8), intent(out) :: gp(ndim, np, nys:nye, nsp)

    integer :: i, j, ii, isp
    real(8) :: dh, sh(-1:1, 2)
    real(8) :: fac1, fac2, txxx, fac1r, fac2r, gam, igam
    real(8) :: bpx, bpy, bpz, epx, epy, epz
    real(8) :: uvm1, uvm2, uvm3, uvm4, uvm5, uvm6
    real(8) :: tmp(1:6, nxs - 1:nxe + 1, nys - 1:nye + 1)

    if (.not. is_init) then
      write (6, *) 'Initialize first by calling particle__init()'
      stop
    end if

    !fields at (i+1/2, j+1/2)
!$OMP PARALLEL DO PRIVATE(i,j)
    do j = nys - 1, nye + 1
    do i = nxs - 1, nxe + 1
      tmp(1, i, j) = 5.0d-1 * (+uf(1, i, j) + uf(1, i, j + 1))
      tmp(2, i, j) = 5.0d-1 * (+uf(2, i, j) + uf(2, i + 1, j))
      tmp(3, i, j) = 2.5d-1 * (+uf(3, i, j) + uf(3, i + 1, j) + uf(3, i, j + 1) + uf(3, i + 1, j + 1))
      tmp(4, i, j) = 5.0d-1 * (+uf(4, i, j) + uf(4, i + 1, j))
      tmp(5, i, j) = 5.0d-1 * (+uf(5, i, j) + uf(5, i, j + 1))
      tmp(6, i, j) = uf(6, i, j)
    end do
    end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO PRIVATE(ii,i,j,isp,sh,dh,gam,igam,fac1,fac2,txxx,fac1r,fac2r, &
!$OMP                     bpx,bpy,bpz,epx,epy,epz,uvm1,uvm2,uvm3,uvm4,uvm5,uvm6)
    do j = nys, nye
    do i = nxs, nxe

      do isp = 1, nsp

        fac1 = q(isp) / r(isp) * 5.0d-1 * delt
        txxx = fac1 * fac1
        fac2 = q(isp) * delt / r(isp)

        do ii = cumcnt(i, j, isp) + 1, cumcnt(i + 1, j, isp)

          !second order shape function
          dh = up(1, ii, j, isp) * d_delx - 5.0d-1 - dble(i)
          sh(-1, 1) = 5.0d-1 * (5.0d-1 - dh) * (5.0d-1 - dh)
          sh(0, 1) = 7.5d-1 - dh * dh
          sh(+1, 1) = 5.0d-1 * (5.0d-1 + dh) * (5.0d-1 + dh)

          dh = up(2, ii, j, isp) * d_delx - 5.0d-1 - dble(j)
          sh(-1, 2) = 5.0d-1 * (5.0d-1 - dh) * (5.0d-1 - dh)
          sh(0, 2) = 7.5d-1 - dh * dh
          sh(+1, 2) = 5.0d-1 * (5.0d-1 + dh) * (5.0d-1 + dh)

          bpx = (tmp(1, i - 1, j - 1) * sh(-1, 1) + tmp(1, i, j - 1) * sh(0, 1) + tmp(1, i + 1, j - 1) * sh(1, 1)) * sh(-1, 2) &
                + (+tmp(1, i - 1, j) * sh(-1, 1) + tmp(1, i, j) * sh(0, 1) + tmp(1, i + 1, j) * sh(+1, 1)) * sh(0, 2) &
                + (+tmp(1, i - 1, j + 1) * sh(-1, 1) + tmp(1, i, j + 1) * sh(0, 1) + tmp(1, i + 1, j + 1) * sh(1, 1)) * sh(1, 2)

          bpy = (tmp(2, i - 1, j - 1) * sh(-1, 1) + tmp(2, i, j - 1) * sh(0, 1) + tmp(2, i + 1, j - 1) * sh(1, 1)) * sh(-1, 2) &
                + (+tmp(2, i - 1, j) * sh(-1, 1) + tmp(2, i, j) * sh(0, 1) + tmp(2, i + 1, j) * sh(+1, 1)) * sh(0, 2) &
                + (+tmp(2, i - 1, j + 1) * sh(-1, 1) + tmp(2, i, j + 1) * sh(0, 1) + tmp(2, i + 1, j + 1) * sh(1, 1)) * sh(1, 2)

          bpz = (tmp(3, i - 1, j - 1) * sh(-1, 1) + tmp(3, i, j - 1) * sh(0, 1) + tmp(3, i + 1, j - 1) * sh(1, 1)) * sh(-1, 2) &
                + (+tmp(3, i - 1, j) * sh(-1, 1) + tmp(3, i, j) * sh(0, 1) + tmp(3, i + 1, j) * sh(+1, 1)) * sh(0, 2) &
                + (+tmp(3, i - 1, j + 1) * sh(-1, 1) + tmp(3, i, j + 1) * sh(0, 1) + tmp(3, i + 1, j + 1) * sh(1, 1)) * sh(1, 2)

          epx = (tmp(4, i - 1, j - 1) * sh(-1, 1) + tmp(4, i, j - 1) * sh(0, 1) + tmp(4, i + 1, j - 1) * sh(1, 1)) * sh(-1, 2) &
                + (+tmp(4, i - 1, j) * sh(-1, 1) + tmp(4, i, j) * sh(0, 1) + tmp(4, i + 1, j) * sh(+1, 1)) * sh(0, 2) &
                + (+tmp(4, i - 1, j + 1) * sh(-1, 1) + tmp(4, i, j + 1) * sh(0, 1) + tmp(4, i + 1, j + 1) * sh(1, 1)) * sh(1, 2)

          epy = (tmp(5, i - 1, j - 1) * sh(-1, 1) + tmp(5, i, j - 1) * sh(0, 1) + tmp(5, i + 1, j - 1) * sh(1, 1)) * sh(-1, 2) &
                + (+tmp(5, i - 1, j) * sh(-1, 1) + tmp(5, i, j) * sh(0, 1) + tmp(5, i + 1, j) * sh(+1, 1)) * sh(0, 2) &
                + (+tmp(5, i - 1, j + 1) * sh(-1, 1) + tmp(5, i, j + 1) * sh(0, 1) + tmp(5, i + 1, j + 1) * sh(1, 1)) * sh(1, 2)

          epz = (tmp(6, i - 1, j - 1) * sh(-1, 1) + tmp(6, i, j - 1) * sh(0, 1) + tmp(6, i + 1, j - 1) * sh(1, 1)) * sh(-1, 2) &
                + (+tmp(6, i - 1, j) * sh(-1, 1) + tmp(6, i, j) * sh(0, 1) + tmp(6, i + 1, j) * sh(+1, 1)) * sh(0, 2) &
                + (+tmp(6, i - 1, j + 1) * sh(-1, 1) + tmp(6, i, j + 1) * sh(0, 1) + tmp(6, i + 1, j + 1) * sh(1, 1)) * sh(1, 2)

          !accel.
          uvm1 = up(3, ii, j, isp) + fac1 * epx
          uvm2 = up(4, ii, j, isp) + fac1 * epy
          uvm3 = up(5, ii, j, isp) + fac1 * epz

          !rotate
          gam = sqrt(c * c + uvm1 * uvm1 + uvm2 * uvm2 + uvm3 * uvm3)
          igam = 1.0d0 / gam
          fac1r = fac1 * igam
          fac2r = fac2 / (gam + txxx * (bpx * bpx + bpy * bpy + bpz * bpz) * igam)

          uvm4 = uvm1 + fac1r * (+uvm2 * bpz - uvm3 * bpy)
          uvm5 = uvm2 + fac1r * (+uvm3 * bpx - uvm1 * bpz)
          uvm6 = uvm3 + fac1r * (+uvm1 * bpy - uvm2 * bpx)

          uvm1 = uvm1 + fac2r * (+uvm5 * bpz - uvm6 * bpy)
          uvm2 = uvm2 + fac2r * (+uvm6 * bpx - uvm4 * bpz)
          uvm3 = uvm3 + fac2r * (+uvm4 * bpy - uvm5 * bpx)

          !accel.
          gp(3, ii, j, isp) = uvm1 + fac1 * epx
          gp(4, ii, j, isp) = uvm2 + fac1 * epy
          gp(5, ii, j, isp) = uvm3 + fac1 * epz

          !move
          gam = 1.0d0 / sqrt(1.0d0 + (gp(3, ii, j, isp) * gp(3, ii, j, isp) &
                                      + gp(4, ii, j, isp) * gp(4, ii, j, isp) &
                                      + gp(5, ii, j, isp) * gp(5, ii, j, isp)) / (c * c))

          gp(1, ii, j, isp) = up(1, ii, j, isp) + gp(3, ii, j, isp) * delt * gam
          gp(2, ii, j, isp) = up(2, ii, j, isp) + gp(4, ii, j, isp) * delt * gam

        end do

      end do

    end do
    end do
!$OMP END PARALLEL DO

    if (ndim == 6) then
!$OMP PARALLEL WORKSHARE
      gp(6, :, :, :) = up(6, :, :, :)
!$OMP END PARALLEL WORKSHARE
    end if

  end subroutine particle__solv

end module particle
