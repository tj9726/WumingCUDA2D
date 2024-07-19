module boundary_periodic
  use mpi
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
  integer :: nerr
  integer, allocatable :: nstat(:)
  real(8), save :: delx, delt, c

contains

  subroutine boundary_periodic__init(ndim_in, np_in, nsp_in, nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in, &
                                     nup_in, ndown_in, mnpi_in, mnpr_in, ncomw_in, nerr_in, nstat_in, &
                                     delx_in, delt_in, c_in)

    integer, intent(in) :: ndim_in, np_in, nsp_in
    integer, intent(in) :: nxgs_in, nxge_in, nygs_in, nyge_in, nys_in, nye_in
    integer, intent(in) :: nup_in, ndown_in, mnpi_in, mnpr_in, ncomw_in, nerr_in, nstat_in(:)
    real(8), intent(in) :: delx_in, delt_in, c_in

    ndim = ndim_in
    np = np_in
    nsp = nsp_in
    nxgs = nxgs_in
    nxge = nxge_in
    nygs = nygs_in
    nyge = nyge_in
    nys = nys_in
    nye = nye_in
    nup = nup_in
    ndown = ndown_in
    mnpi = mnpi_in
    mnpr = mnpr_in
    ncomw = ncomw_in
    nerr = nerr_in
    allocate (nstat(size(nstat_in)))
    nstat = nstat_in
    delx = delx_in
    delt = delt_in
    c = c_in

    is_init = .true.

  end subroutine boundary_periodic__init

  subroutine boundary_periodic__particle_x(up, np2)

    use, intrinsic :: ieee_arithmetic

    integer, intent(in) :: np2(nys:nye, nsp)
    real(8), intent(inout) :: up(ndim, np, nys:nye, nsp)
    integer :: j, ii, isp, ipos

    if (.not. is_init) then
      write (6, *) 'Initialize first by calling boundary_periodic__init()'
      stop
    end if

    call ieee_set_rounding_mode(ieee_down)

    do isp = 1, nsp

!$OMP PARALLEL DO PRIVATE(ii,j,ipos)
      do j = nys, nye
        do ii = 1, np2(j, isp)

          ipos = int(up(1, ii, j, isp) / delx)

          if (ipos < nxgs) then
            up(1, ii, j, isp) = up(1, ii, j, isp) + (nxge - nxgs + 1) * delx
          else if (ipos >= nxge + 1) then
            up(1, ii, j, isp) = up(1, ii, j, isp) - (nxge - nxgs + 1) * delx
          end if

        end do
      end do
!$OMP END PARALLEL DO

    end do

  end subroutine boundary_periodic__particle_x

  subroutine boundary_periodic__particle_y(up, np2)

    use, intrinsic :: ieee_arithmetic
!$  use omp_lib

    integer, intent(inout) :: np2(nys:nye, nsp)
    real(8), intent(inout) :: up(ndim, np, nys:nye, nsp)
    logical, save :: lflag = .true.
!$  integer(omp_lock_kind) :: lck(nys - 1:nye + 1)
    integer :: j, ii, iii, isp, jpos, idim
    integer :: cnt(nys - 1:nye + 1), cnt2(nys:nye), cnt_tmp
    integer, save, allocatable :: flag(:, :)
    real(8), save, allocatable :: bff_ptcl(:, :)

    if (.not. is_init) then
      write (6, *) 'Initialize first by calling boundary_periodic__init()'
      stop
    end if

    if (lflag) then
      allocate (flag(np, nys:nye))
      allocate (bff_ptcl(ndim * np, nys - 1:nye + 1))
      lflag = .false.
    end if

    call ieee_set_rounding_mode(ieee_down)

!$OMP PARALLEL DO PRIVATE(j)
!$  do j = nys - 1, nye + 1
!$    call omp_init_lock(lck(j))
!$  end do
!$OMP END PARALLEL DO

    do isp = 1, nsp

!$OMP PARALLEL

!$OMP WORKSHARE
      cnt(nys - 1:nye + 1) = 0
!$OMP END WORKSHARE
!$OMP WORKSHARE
      cnt2(nys:nye) = 0
!$OMP END WORKSHARE

!$OMP DO PRIVATE(ii,j,idim,jpos)
      do j = nys, nye
        do ii = 1, np2(j, isp)

          jpos = int(up(2, ii, j, isp) / delx)

          if (jpos /= j) then
            if (jpos <= nygs - 1) then
              up(2, ii, j, isp) = up(2, ii, j, isp) + (nyge - nygs + 1) * delx
            else if (jpos >= nyge + 1) then
              up(2, ii, j, isp) = up(2, ii, j, isp) - (nyge - nygs + 1) * delx
            end if

!$          call omp_set_lock(lck(jpos))
            do idim = 1, ndim
              bff_ptcl(idim + ndim * cnt(jpos), jpos) = up(idim, ii, j, isp)
            end do
            cnt(jpos) = cnt(jpos) + 1
!$          call omp_unset_lock(lck(jpos))

            cnt2(j) = cnt2(j) + 1
            flag(cnt2(j), j) = ii
          end if

        end do
      end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

      !transfer to rank-1
      call MPI_SENDRECV(cnt(nys - 1), 1, mnpi, ndown, 100, &
                        cnt_tmp, 1, mnpi, nup, 100, &
                        ncomw, nstat, nerr)
      call MPI_SENDRECV(bff_ptcl(1, nys - 1), ndim * cnt(nys - 1), mnpr, ndown, 101, &
                        bff_ptcl(ndim * cnt(nye) + 1, nye), ndim * cnt_tmp, mnpr, nup, 101, &
                        ncomw, nstat, nerr)
      cnt(nye) = cnt(nye) + cnt_tmp

      !transfer to rank+1
      call MPI_SENDRECV(cnt(nye + 1), 1, mnpi, nup, 200, &
                        cnt_tmp, 1, mnpi, ndown, 200, &
                        ncomw, nstat, nerr)
      call MPI_SENDRECV(bff_ptcl(1, nye + 1), ndim * cnt(nye + 1), mnpr, nup, 201, &
                        bff_ptcl(ndim * cnt(nys) + 1, nys), ndim * cnt_tmp, mnpr, ndown, 201, &
                        ncomw, nstat, nerr)
      cnt(nys) = cnt(nys) + cnt_tmp

!$OMP PARALLEL

!$OMP DO PRIVATE(iii,ii,j,idim,cnt_tmp)
      do j = nys, nye
        iii = 0
        cnt_tmp = cnt2(j)
        loop1: do ii = 1, cnt2(j)
          if (cnt(j) == 0) then
            if (np2(j, isp) < flag(ii, j)) exit loop1
            do while (np2(j, isp) == flag(cnt_tmp, j))
              np2(j, isp) = np2(j, isp) - 1
              if (np2(j, isp) < flag(ii, j)) exit loop1
              cnt_tmp = cnt_tmp - 1
            end do
            do idim = 1, ndim
              up(idim, flag(ii, j), j, isp) = up(idim, np2(j, isp), j, isp)
            end do
            np2(j, isp) = np2(j, isp) - 1
          else
            do idim = 1, ndim
              up(idim, flag(ii, j), j, isp) = bff_ptcl(idim + ndim * iii, j)
            end do
            iii = iii + 1
            cnt(j) = cnt(j) - 1
          end if
        end do loop1

        if (cnt(j) > 0) then
          do ii = 1, cnt(j)
            do idim = 1, ndim
              up(idim, np2(j, isp) + ii, j, isp) = bff_ptcl(ndim * iii + idim + ndim * (ii - 1), j)
            end do
          end do
        end if
      end do
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(j)
      do j = nys, nye
        np2(j, isp) = np2(j, isp) + cnt(j)
        if (np2(j, isp) > np) then
          write (*, *) "memory over (np2 > np)", np, np2(j, isp), j, isp
          stop
        end if
      end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    end do

!$OMP PARALLEL DO PRIVATE(j)
!$  do j = nys - 1, nye + 1
!$    call omp_destroy_lock(lck(j))
!$  end do
!$OMP END PARALLEL DO

  end subroutine boundary_periodic__particle_y

  subroutine boundary_periodic__dfield(df, nxs, nxe, nys, nye, nxgs, nxge)

    integer, intent(in) :: nxs, nxe, nys, nye, nxgs, nxge
    real(8), intent(inout) :: df(6, nxgs - 2:nxge + 2, nys - 2:nye + 2)
    integer :: i, ii
    real(8) :: bff_snd(12 * (nxe - nxs + 1)), bff_rcv(12 * (nxe - nxs + 1))

    if (.not. is_init) then
      write (6, *) 'Initialize first by calling boundary_periodic__init()'
      stop
    end if

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = 12 * (i - nxs)
      bff_snd(ii + 1) = df(1, i, nys)
      bff_snd(ii + 2) = df(2, i, nys)
      bff_snd(ii + 3) = df(3, i, nys)
      bff_snd(ii + 4) = df(4, i, nys)
      bff_snd(ii + 5) = df(5, i, nys)
      bff_snd(ii + 6) = df(6, i, nys)
      bff_snd(ii + 7) = df(1, i, nys + 1)
      bff_snd(ii + 8) = df(2, i, nys + 1)
      bff_snd(ii + 9) = df(3, i, nys + 1)
      bff_snd(ii + 10) = df(4, i, nys + 1)
      bff_snd(ii + 11) = df(5, i, nys + 1)
      bff_snd(ii + 12) = df(6, i, nys + 1)
    end do
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1), 12 * (nxe - nxs + 1), mnpr, ndown, 110, &
                      bff_rcv(1), 12 * (nxe - nxs + 1), mnpr, nup, 110, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = 12 * (i - nxs)
      df(1, i, nye + 1) = bff_rcv(ii + 1)
      df(2, i, nye + 1) = bff_rcv(ii + 2)
      df(3, i, nye + 1) = bff_rcv(ii + 3)
      df(4, i, nye + 1) = bff_rcv(ii + 4)
      df(5, i, nye + 1) = bff_rcv(ii + 5)
      df(6, i, nye + 1) = bff_rcv(ii + 6)
      df(1, i, nye + 2) = bff_rcv(ii + 7)
      df(2, i, nye + 2) = bff_rcv(ii + 8)
      df(3, i, nye + 2) = bff_rcv(ii + 9)
      df(4, i, nye + 2) = bff_rcv(ii + 10)
      df(5, i, nye + 2) = bff_rcv(ii + 11)
      df(6, i, nye + 2) = bff_rcv(ii + 12)
    end do
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = 12 * (i - nxs)
      bff_snd(ii + 1) = df(1, i, nye - 1)
      bff_snd(ii + 2) = df(2, i, nye - 1)
      bff_snd(ii + 3) = df(3, i, nye - 1)
      bff_snd(ii + 4) = df(4, i, nye - 1)
      bff_snd(ii + 5) = df(5, i, nye - 1)
      bff_snd(ii + 6) = df(6, i, nye - 1)
      bff_snd(ii + 7) = df(1, i, nye)
      bff_snd(ii + 8) = df(2, i, nye)
      bff_snd(ii + 9) = df(3, i, nye)
      bff_snd(ii + 10) = df(4, i, nye)
      bff_snd(ii + 11) = df(5, i, nye)
      bff_snd(ii + 12) = df(6, i, nye)
    end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1), 12 * (nxe - nxs + 1), mnpr, nup, 100, &
                      bff_rcv(1), 12 * (nxe - nxs + 1), mnpr, ndown, 100, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = 12 * (i - nxs)
      df(1, i, nys - 2) = bff_rcv(ii + 1)
      df(2, i, nys - 2) = bff_rcv(ii + 2)
      df(3, i, nys - 2) = bff_rcv(ii + 3)
      df(4, i, nys - 2) = bff_rcv(ii + 4)
      df(5, i, nys - 2) = bff_rcv(ii + 5)
      df(6, i, nys - 2) = bff_rcv(ii + 6)
      df(1, i, nys - 1) = bff_rcv(ii + 7)
      df(2, i, nys - 1) = bff_rcv(ii + 8)
      df(3, i, nys - 1) = bff_rcv(ii + 9)
      df(4, i, nys - 1) = bff_rcv(ii + 10)
      df(5, i, nys - 1) = bff_rcv(ii + 11)
      df(6, i, nys - 1) = bff_rcv(ii + 12)
    end do
!$OMP END PARALLEL DO

!$OMP WORKSHARE
    df(1:6, nxs - 2, nys - 2:nye + 2) = df(1:6, nxe - 1, nys - 2:nye + 2)
    df(1:6, nxs - 1, nys - 2:nye + 2) = df(1:6, nxe, nys - 2:nye + 2)
    df(1:6, nxe + 1, nys - 2:nye + 2) = df(1:6, nxs, nys - 2:nye + 2)
    df(1:6, nxe + 2, nys - 2:nye + 2) = df(1:6, nxs + 1, nys - 2:nye + 2)
!$OMP END WORKSHARE

  end subroutine boundary_periodic__dfield

  subroutine boundary_periodic__curre(uj, nxs, nxe, nys, nye, nxgs, nxge)

    integer, intent(in) :: nxs, nxe, nys, nye, nxgs, nxge
    real(8), intent(inout) :: uj(3, nxgs - 2:nxge + 2, nys - 2:nye + 2)
    integer :: i, ii
    real(8) :: bff_rcv(6 * (nxe - nxs + 4 + 1)), bff_snd(6 * (nxe - nxs + 4 + 1))

    if (.not. is_init) then
      write (6, *) 'Initialize first by calling boundary_periodic__init()'
      stop
    end if

    !send to rank-1
!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      bff_snd(ii + 1) = uj(1, i, nys - 2)
      bff_snd(ii + 2) = uj(2, i, nys - 2)
      bff_snd(ii + 3) = uj(3, i, nys - 2)
      bff_snd(ii + 4) = uj(1, i, nys - 1)
      bff_snd(ii + 5) = uj(2, i, nys - 1)
      bff_snd(ii + 6) = uj(3, i, nys - 1)
    end do
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1), 6 * (nxe - nxs + 4 + 1), mnpr, ndown, 110, &
                      bff_rcv(1), 6 * (nxe - nxs + 4 + 1), mnpr, nup, 110, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      uj(1, i, nye - 1) = uj(1, i, nye - 1) + bff_rcv(ii + 1)
      uj(2, i, nye - 1) = uj(2, i, nye - 1) + bff_rcv(ii + 2)
      uj(3, i, nye - 1) = uj(3, i, nye - 1) + bff_rcv(ii + 3)
      uj(1, i, nye) = uj(1, i, nye) + bff_rcv(ii + 4)
      uj(2, i, nye) = uj(2, i, nye) + bff_rcv(ii + 5)
      uj(3, i, nye) = uj(3, i, nye) + bff_rcv(ii + 6)
    end do
!$OMP END DO NOWAIT

    !send to rank+1
!$OMP DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      bff_snd(ii + 1) = uj(1, i, nye + 1)
      bff_snd(ii + 2) = uj(2, i, nye + 1)
      bff_snd(ii + 3) = uj(3, i, nye + 1)
      bff_snd(ii + 4) = uj(1, i, nye + 2)
      bff_snd(ii + 5) = uj(2, i, nye + 2)
      bff_snd(ii + 6) = uj(3, i, nye + 2)
    end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1), 6 * (nxe - nxs + 4 + 1), mnpr, nup, 120, &
                      bff_rcv(1), 6 * (nxe - nxs + 4 + 1), mnpr, ndown, 120, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      uj(1, i, nys) = uj(1, i, nys) + bff_rcv(ii + 1)
      uj(2, i, nys) = uj(2, i, nys) + bff_rcv(ii + 2)
      uj(3, i, nys) = uj(3, i, nys) + bff_rcv(ii + 3)
      uj(1, i, nys + 1) = uj(1, i, nys + 1) + bff_rcv(ii + 4)
      uj(2, i, nys + 1) = uj(2, i, nys + 1) + bff_rcv(ii + 5)
      uj(3, i, nys + 1) = uj(3, i, nys + 1) + bff_rcv(ii + 6)
    end do
!$OMP END PARALLEL DO

!#####    !Update of nori-shiro   #####

    !send to rank-1
!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      bff_snd(ii + 1) = uj(1, i, nys)
      bff_snd(ii + 2) = uj(2, i, nys)
      bff_snd(ii + 3) = uj(3, i, nys)
      bff_snd(ii + 4) = uj(1, i, nys + 1)
      bff_snd(ii + 5) = uj(2, i, nys + 1)
      bff_snd(ii + 6) = uj(3, i, nys + 1)
    end do
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1), 6 * (nxe - nxs + 4 + 1), mnpr, ndown, 130, &
                      bff_rcv(1), 6 * (nxe - nxs + 4 + 1), mnpr, nup, 130, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      uj(1, i, nye + 1) = bff_rcv(ii + 1)
      uj(2, i, nye + 1) = bff_rcv(ii + 2)
      uj(3, i, nye + 1) = bff_rcv(ii + 3)
      uj(1, i, nye + 2) = bff_rcv(ii + 4)
      uj(2, i, nye + 2) = bff_rcv(ii + 5)
      uj(3, i, nye + 2) = bff_rcv(ii + 6)
    end do
!$OMP END DO NOWAIT

    !send to rank+1
!$OMP DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      bff_snd(ii + 1) = uj(1, i, nye - 1)
      bff_snd(ii + 2) = uj(2, i, nye - 1)
      bff_snd(ii + 3) = uj(3, i, nye - 1)
      bff_snd(ii + 4) = uj(1, i, nye)
      bff_snd(ii + 5) = uj(2, i, nye)
      bff_snd(ii + 6) = uj(3, i, nye)
    end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1), 6 * (nxe - nxs + 4 + 1), mnpr, nup, 140, &
                      bff_rcv(1), 6 * (nxe - nxs + 4 + 1), mnpr, ndown, 140, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs - 2, nxe + 2
      ii = 6 * (i - (nxs - 2))
      uj(1, i, nys - 2) = bff_rcv(ii + 1)
      uj(2, i, nys - 2) = bff_rcv(ii + 2)
      uj(3, i, nys - 2) = bff_rcv(ii + 3)
      uj(1, i, nys - 1) = bff_rcv(ii + 4)
      uj(2, i, nys - 1) = bff_rcv(ii + 5)
      uj(3, i, nys - 1) = bff_rcv(ii + 6)
    end do
!$OMP END PARALLEL DO

!$OMP WORKSHARE
    uj(1:3, nxe - 1, nys - 2:nye + 2) = uj(1:3, nxe - 1, nys - 2:nye + 2) + uj(1:3, nxs - 2, nys - 2:nye + 2)
    uj(1:3, nxe, nys - 2:nye + 2) = uj(1:3, nxe, nys - 2:nye + 2) + uj(1:3, nxs - 1, nys - 2:nye + 2)
    uj(1:3, nxs, nys - 2:nye + 2) = uj(1:3, nxs, nys - 2:nye + 2) + uj(1:3, nxe + 1, nys - 2:nye + 2)
    uj(1:3, nxs + 1, nys - 2:nye + 2) = uj(1:3, nxs + 1, nys - 2:nye + 2) + uj(1:3, nxe + 2, nys - 2:nye + 2)
!$OMP END WORKSHARE
!$OMP WORKSHARE
    uj(1:3, nxs - 2, nys - 2:nye + 2) = uj(1:3, nxe - 1, nys - 2:nye + 2)
    uj(1:3, nxs - 1, nys - 2:nye + 2) = uj(1:3, nxe, nys - 2:nye + 2)
    uj(1:3, nxe + 1, nys - 2:nye + 2) = uj(1:3, nxs, nys - 2:nye + 2)
    uj(1:3, nxe + 2, nys - 2:nye + 2) = uj(1:3, nxs + 1, nys - 2:nye + 2)
!$OMP END WORKSHARE

  end subroutine boundary_periodic__curre

  subroutine boundary_periodic__phi(phi, nxs, nxe, nys, nye, l)

    integer, intent(in) :: nxs, nxe, nys, nye, l
    real(8), intent(inout) :: phi(nxs - 1:nxe + 1, nys - 1:nye + 1)
    integer :: i, ii
    real(8) :: bff_snd(nxe - nxs + 1), bff_rcv(nxe - nxs + 1)

    if (.not. is_init) then
      write (6, *) 'Initialize first by calling boundary_periodic__init()'
      stop
    end if

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = i - nxs
      bff_snd(ii + 1) = phi(i, nys)
    end do
!$OMP END PARALLEL DO

    call MPI_SENDRECV(bff_snd(1), nxe - nxs + 1, mnpr, ndown, 110, &
                      bff_rcv(1), nxe - nxs + 1, mnpr, nup, 110, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL

!$OMP DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = i - nxs
      phi(i, nye + 1) = bff_rcv(ii + 1)
    end do
!$OMP END DO NOWAIT

!$OMP DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = i - nxs
      bff_snd(ii + 1) = phi(i, nye)
    end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    call MPI_SENDRECV(bff_snd(1), nxe - nxs + 1, mnpr, nup, 100, &
                      bff_rcv(1), nxe - nxs + 1, mnpr, ndown, 100, &
                      ncomw, nstat, nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
    do i = nxs, nxe
      ii = i - nxs
      phi(i, nys - 1) = bff_rcv(ii + 1)
    end do
!$OMP END PARALLEL DO

!$OMP WORKSHARE
    phi(nxs - 1, nys - 1:nye + 1) = phi(nxe, nys - 1:nye + 1)
    phi(nxe + 1, nys - 1:nye + 1) = phi(nxs, nys - 1:nye + 1)
!$OMP END WORKSHARE

  end subroutine boundary_periodic__phi

  subroutine boundary_periodic__mom(nvector, ttensor)

    real(8), intent(inout) :: nvector(4, nxgs - 1:nxge + 1, nys - 1:nye + 1, nsp), &
                              ttensor(10, nxgs - 1:nxge + 1, nys - 1:nye + 1, nsp)

    integer :: i, ii, k, isp
    real(8) :: bff_rcv(10 * (nxge - nxgs + 3)), bff_snd(10 * (nxge - nxgs + 3))

!$OMP PARALLEL WORKSHARE
    nvector(1:4, nxgs, nys - 1:nye + 1, 1:nsp) = nvector(1:4, nxgs, nys - 1:nye + 1, 1:nsp) &
                                                 + nvector(1:4, nxge + 1, nys - 1:nye + 1, 1:nsp)
    nvector(1:4, nxge, nys - 1:nye + 1, 1:nsp) = nvector(1:4, nxge, nys - 1:nye + 1, 1:nsp) &
                                                 + nvector(1:4, nxgs - 1, nys - 1:nye + 1, 1:nsp)
    ttensor(1:10, nxgs, nys - 1:nye + 1, 1:nsp) = ttensor(1:10, nxgs, nys - 1:nye + 1, 1:nsp) &
                                                  + ttensor(1:10, nxge + 1, nys - 1:nye + 1, 1:nsp)
    ttensor(1:10, nxge, nys - 1:nye + 1, 1:nsp) = ttensor(1:10, nxge, nys - 1:nye + 1, 1:nsp) &
                                                  + ttensor(1:10, nxgs - 1, nys - 1:nye + 1, 1:nsp)
!$OMP END PARALLEL WORKSHARE

    !nvector
    do isp = 1, nsp
      !send to rank-1
!$OMP PARALLEL DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 4 * (i - (nxgs - 1))
        do k = 1, 4
          bff_snd(ii + k) = nvector(k, i, nys - 1, isp)
        end do
      end do
!$OMP END PARALLEL DO

      call MPI_SENDRECV(bff_snd(1), 4 * (nxge - nxgs + 3), mnpr, ndown, 200, &
                        bff_rcv(1), 4 * (nxge - nxgs + 3), mnpr, nup, 200, &
                        ncomw, nstat, nerr)

!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 4 * (i - (nxgs - 1))
        do k = 1, 4
          nvector(k, i, nye, isp) = nvector(k, i, nye, isp) + bff_rcv(ii + k)
        end do
      end do
!$OMP END DO NOWAIT

      !send to rank+1
!$OMP DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 4 * (i - (nxgs - 1))
        do k = 1, 4
          bff_snd(ii + k) = nvector(k, i, nye + 1, isp)
        end do
      end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      call MPI_SENDRECV(bff_snd(1), 4 * (nxge - nxgs + 3), mnpr, nup, 201, &
                        bff_rcv(1), 4 * (nxge - nxgs + 3), mnpr, ndown, 201, &
                        ncomw, nstat, nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 4 * (i - (nxgs - 1))
        do k = 1, 4
          nvector(k, i, nys, isp) = nvector(k, i, nys, isp) + bff_rcv(ii + k)
        end do
      end do
!$OMP END PARALLEL DO
    end do

    !ttensor
    do isp = 1, nsp
      !send to rank-1
!$OMP PARALLEL DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 10 * (i - (nxgs - 1))
        do k = 1, 10
          bff_snd(ii + k) = ttensor(k, i, nys - 1, isp)
        end do
      end do
!$OMP END PARALLEL DO

      call MPI_SENDRECV(bff_snd(1), 10 * (nxge - nxgs + 3), mnpr, ndown, 200, &
                        bff_rcv(1), 10 * (nxge - nxgs + 3), mnpr, nup, 200, &
                        ncomw, nstat, nerr)

!$OMP PARALLEL
!$OMP DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 10 * (i - (nxgs - 1))
        do k = 1, 10
          ttensor(k, i, nye, isp) = ttensor(k, i, nye, isp) + bff_rcv(ii + k)
        end do
      end do
!$OMP END DO NOWAIT

      !send to rank+1
!$OMP DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 10 * (i - (nxgs - 1))
        do k = 1, 10
          bff_snd(ii + k) = ttensor(k, i, nye + 1, isp)
        end do
      end do
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      call MPI_SENDRECV(bff_snd(1), 10 * (nxge - nxgs + 3), mnpr, nup, 201, &
                        bff_rcv(1), 10 * (nxge - nxgs + 3), mnpr, ndown, 201, &
                        ncomw, nstat, nerr)

!$OMP PARALLEL DO PRIVATE(i,ii)
      do i = nxgs - 1, nxge + 1
        ii = 10 * (i - (nxgs - 1))
        do k = 1, 10
          ttensor(k, i, nys, isp) = ttensor(k, i, nys, isp) + bff_rcv(ii + k)
        end do
      end do
!$OMP END PARALLEL DO
    end do

  end subroutine boundary_periodic__mom

end module boundary_periodic
