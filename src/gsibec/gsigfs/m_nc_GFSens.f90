module m_nc_GFSens
!$$$ module documentation block
!                .      .    .                                       .
! module:   m_nc_GFSens
!   prgmmr: todling          org: np22                date: 2024-01-01
!
! abstract: NetCDF4 I/O routines for GFS gaussian grid ensemble files.
!           Analogous to m_nc_GEOSens.f90 but adapted for GFS conventions.
!
! program history log:
!   2024-01-01  initial version following m_nc_GEOSens.f90 structure
!
! subroutines included:
!   sub nc_GFSens_dims       - read file dimensions
!   sub nc_GFSens_read       - read GFS ensemble member from NetCDF4 file
!   sub nc_GFSens_vars_set   - initialize nc_GFSens_vars type from CV names
!   sub nc_GFSens_vars_init  - allocate arrays in nc_GFSens_vars
!   sub nc_GFSens_vars_final - deallocate arrays in nc_GFSens_vars
!   sub nc_GFSens_summary    - print summary of nc_GFSens_vars
!   sub nc_GFSens_gfs2gsi    - convert from GFS to GSI units/orientation
!
! attributes:
!   language: Fortran 90 and/or above
!
!$$$

use netcdf
use mpeu_util, only: getindex
implicit none
private

public :: nc_GFSens_vars_set
public :: nc_GFSens_vars_init
public :: nc_GFSens_vars_final
public :: nc_GFSens_vars_comp
public :: nc_GFSens_vars
public :: nc_GFSens_dims
public :: nc_GFSens_read
public :: nc_GFSens_summary
public :: nc_GFSens_gfs2gsi
public :: nc_GFSens_getpointer

! This type stores GFS ensemble data for one member.
! Arrays are stored in GSI orientation (nlat, nlon, [nsig,] nvar)
! when gsiset=.true., else (nlon, nlat, [nsig,] nvar).
type nc_GFSens_vars
   logical :: initialized = .false.
   integer :: nlon = 0, nlat = 0, nsig = 0
   logical :: gsiset = .false.
   real(4), pointer, dimension(:,:,:,:) :: ptr3d => null()  ! (nlat,nlon,nsig,nv3d) if gsiset
   real(4), pointer, dimension(:,:,:)   :: ptr2d => null()  ! (nlat,nlon,nv2d) if gsiset
   integer :: nv2d = -1
   integer :: nv3d = -1
   character(len=8), allocatable :: gsi_vnames2d(:)  ! GSI control variable names (2d)
   character(len=8), allocatable :: gsi_vnames3d(:)  ! GSI control variable names (3d)
end type nc_GFSens_vars

character(len=*), parameter :: myname = 'm_nc_GFSens'

! Unit conversion: surface pressure Pa -> centibars (cb = kPa)
! 1 cb = 1 kPa = 1000 Pa
real, parameter :: Pa_to_cb = 1.0e-3  ! Pa to centibars

interface nc_GFSens_dims;      module procedure read_dims_;         end interface
interface nc_GFSens_read;      module procedure read_GFSens_;       end interface
interface nc_GFSens_vars_set;  module procedure set_vars_;          end interface
interface nc_GFSens_vars_init; module procedure init_GFSens_vars_;  end interface
interface nc_GFSens_vars_final;module procedure final_GFSens_vars_; end interface
interface nc_GFSens_vars_comp; module procedure comp_GFSens_vars_;  end interface
interface nc_GFSens_summary;   module procedure summary_;           end interface
interface nc_GFSens_gfs2gsi;   module procedure gfs2gsi_;           end interface
interface nc_GFSens_getpointer
   module procedure get_pointer_2d_
   module procedure get_pointer_3d_
end interface

contains

!---------------------------------------------------------------------------
! Return the GFS file variable name corresponding to a GSI control variable name.
! GFS NetCDF4 gaussian grid ensemble files use standard NCEP variable names.
function gfs_varname_(gsiname) result(gfsname)
   character(len=*), intent(in) :: gsiname
   character(len=32) :: gfsname
   select case(trim(adjustl(gsiname)))
   case('sf')   ; gfsname = 'ugrd'    ! zonal wind -> stream function slot (uv_hyb_ens)
   case('vp')   ; gfsname = 'vgrd'    ! meridional wind -> vel. potential slot (uv_hyb_ens)
   case('u')    ; gfsname = 'ugrd'
   case('v')    ; gfsname = 'vgrd'
   case('t')    ; gfsname = 'tmp'     ! temperature
   case('tv')   ; gfsname = 'tmp'     ! virtual temperature (file has T, GSI converts)
   case('q')    ; gfsname = 'spfh'    ! specific humidity
   case('oz')   ; gfsname = 'o3mr'    ! ozone mixing ratio
   case('cw')   ; gfsname = 'clwmr'  ! cloud liquid water
   case('ql')   ; gfsname = 'clwmr'  ! cloud liquid water
   case('qi')   ; gfsname = 'icmr'   ! cloud ice
   case('qr')   ; gfsname = 'rwmr'   ! rain water
   case('qs')   ; gfsname = 'snmr'   ! snow
   case('qg')   ; gfsname = 'grle'   ! graupel
   case('ps')   ; gfsname = 'pressfc' ! surface pressure
   case default ; gfsname = trim(gsiname)
   end select
end function gfs_varname_

!---------------------------------------------------------------------------
subroutine set_vars_(cvars2d, cvars3d, bvars)
   implicit none
   character(*), intent(in)  :: cvars2d(:)  ! GSI 2d control variable names
   character(*), intent(in)  :: cvars3d(:)  ! GSI 3d control variable names
   type(nc_GFSens_vars), intent(inout) :: bvars

   bvars%nv2d = size(cvars2d)
   bvars%nv3d = size(cvars3d)
   if (bvars%nv2d > 0) then
      if (.not. allocated(bvars%gsi_vnames2d)) &
           allocate(bvars%gsi_vnames2d(bvars%nv2d))
      bvars%gsi_vnames2d = cvars2d
   endif
   if (bvars%nv3d > 0) then
      if (.not. allocated(bvars%gsi_vnames3d)) &
           allocate(bvars%gsi_vnames3d(bvars%nv3d))
      bvars%gsi_vnames3d = cvars3d
   endif
end subroutine set_vars_

!---------------------------------------------------------------------------
subroutine init_GFSens_vars_(bvars, nlon, nlat, nsig, gsi)
   implicit none
   type(nc_GFSens_vars), intent(inout) :: bvars
   integer, intent(in) :: nlon, nlat, nsig
   logical, intent(in), optional :: gsi

   bvars%nlon = nlon
   bvars%nlat = nlat
   bvars%nsig = nsig
   bvars%initialized = .true.

   if (present(gsi)) then
      bvars%gsiset = gsi
   else
      bvars%gsiset = .false.
   endif

   if (bvars%nv2d > 0) then
      if (bvars%gsiset) then
         allocate(bvars%ptr2d(nlat, nlon, bvars%nv2d))
      else
         allocate(bvars%ptr2d(nlon, nlat, bvars%nv2d))
      endif
   endif
   if (bvars%nv3d > 0) then
      if (bvars%gsiset) then
         allocate(bvars%ptr3d(nlat, nlon, nsig, bvars%nv3d))
      else
         allocate(bvars%ptr3d(nlon, nlat, nsig, bvars%nv3d))
      endif
   endif
end subroutine init_GFSens_vars_

!---------------------------------------------------------------------------
subroutine final_GFSens_vars_(bvars)
   implicit none
   type(nc_GFSens_vars), intent(inout) :: bvars

   ! Deallocate data arrays and reset initialized flag.
   ! NOTE: gsi_vnames2d/gsi_vnames3d are preserved for reuse on next call.
   bvars%initialized = .false.
   if (bvars%nv2d > 0 .and. associated(bvars%ptr2d)) then
      deallocate(bvars%ptr2d)
      nullify(bvars%ptr2d)
   endif
   if (bvars%nv3d > 0 .and. associated(bvars%ptr3d)) then
      deallocate(bvars%ptr3d)
      nullify(bvars%ptr3d)
   endif
end subroutine final_GFSens_vars_

!---------------------------------------------------------------------------
subroutine comp_GFSens_vars_(avars, bvars, verbose)
   implicit none
   type(nc_GFSens_vars), intent(in) :: avars, bvars
   logical, intent(in), optional :: verbose
   logical :: verbose_
   verbose_ = .false.
   if (present(verbose)) verbose_ = verbose
   if (avars%nlon /= bvars%nlon .or. avars%nlat /= bvars%nlat .or. &
       avars%nsig /= bvars%nsig) then
      if (verbose_) print *, myname, ': fields are inconsistent'
   else
      if (verbose_) print *, myname, ': fields match'
   endif
end subroutine comp_GFSens_vars_

!---------------------------------------------------------------------------
subroutine summary_(bvars)
   implicit none
   type(nc_GFSens_vars), intent(in) :: bvars
   integer :: nv
   print *, myname, ': nlon,nlat,nsig = ', bvars%nlon, bvars%nlat, bvars%nsig
   print *, myname, ': nv2d,nv3d = ', bvars%nv2d, bvars%nv3d
   print *, myname, ': gsiset = ', bvars%gsiset
   if (allocated(bvars%gsi_vnames2d)) then
      do nv = 1, bvars%nv2d
         print *, myname, ': 2d var ', nv, ' = ', trim(bvars%gsi_vnames2d(nv))
      enddo
   endif
   if (allocated(bvars%gsi_vnames3d)) then
      do nv = 1, bvars%nv3d
         print *, myname, ': 3d var ', nv, ' = ', trim(bvars%gsi_vnames3d(nv))
      enddo
   endif
end subroutine summary_

!---------------------------------------------------------------------------
subroutine read_dims_(fname, nlat, nlon, nlev, rc, myid, root)
   implicit none
   character(len=*), intent(in)  :: fname
   integer, intent(out) :: nlat, nlon, nlev, rc
   integer, intent(in), optional :: myid, root

   integer :: ncid, varid, ier
   integer :: mype_, root_
   character(len=*), parameter :: myname_ = myname//'::read_dims_'

   rc = 0; mype_ = 0; root_ = 0
   if (present(myid) .and. present(root)) then
      mype_ = myid
      root_ = root
   endif

   call check_(nf90_open(fname, NF90_NOWRITE, ncid), rc, mype_, root_)
   if (rc /= 0) return

   ! Try various dimension naming conventions
   ! Attempt 'lon' first (GEOS convention), then 'grid_xt' (FV3GFS convention)
   nlon = 0
   ier = nf90_inq_dimid(ncid, 'lon', varid)
   if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlon), rc, mype_, root_)
   if (nlon == 0) then
      ier = nf90_inq_dimid(ncid, 'grid_xt', varid)
      if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlon), rc, mype_, root_)
   endif
   if (nlon == 0) then
      ier = nf90_inq_dimid(ncid, 'longitude', varid)
      if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlon), rc, mype_, root_)
   endif

   nlat = 0
   ier = nf90_inq_dimid(ncid, 'lat', varid)
   if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlat), rc, mype_, root_)
   if (nlat == 0) then
      ier = nf90_inq_dimid(ncid, 'grid_yt', varid)
      if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlat), rc, mype_, root_)
   endif
   if (nlat == 0) then
      ier = nf90_inq_dimid(ncid, 'latitude', varid)
      if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlat), rc, mype_, root_)
   endif

   nlev = 0
   ier = nf90_inq_dimid(ncid, 'lev', varid)
   if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlev), rc, mype_, root_)
   if (nlev == 0) then
      ier = nf90_inq_dimid(ncid, 'pfull', varid)
      if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlev), rc, mype_, root_)
   endif
   if (nlev == 0) then
      ier = nf90_inq_dimid(ncid, 'level', varid)
      if (ier == NF90_NOERR) call check_(nf90_inquire_dimension(ncid, varid, len=nlev), rc, mype_, root_)
   endif

   if (nlon == 0 .or. nlat == 0 .or. nlev == 0) then
      if (mype_ == root_) then
         print *, myname_, ': could not determine dimensions, nlon,nlat,nlev = ', nlon, nlat, nlev
      endif
      rc = 99
   endif

   call check_(nf90_close(ncid), rc, mype_, root_)
end subroutine read_dims_

!---------------------------------------------------------------------------
subroutine read_GFSens_(fname, bvars, rc, myid, root, gsiset)
   implicit none
   character(len=*), intent(in)    :: fname
   type(nc_GFSens_vars), intent(inout) :: bvars
   integer, intent(out) :: rc
   integer, intent(in), optional :: myid, root
   logical, intent(in), optional :: gsiset

   integer :: ncid, varid, ier
   integer :: kk, nv, nlat, nlon, nlev
   integer :: nlat_, nlon_, nlev_
   integer :: mype_, root_
   real(4), allocatable :: data3d(:,:,:)
   real(4), allocatable :: data2d(:,:)
   logical :: gsi_, verbose, init_
   character(len=32) :: fvname
   character(len=*), parameter :: myname_ = myname//'::read_GFSens_'

   rc = 0; mype_ = 0; root_ = 0
   verbose = .true.
   if (present(myid) .and. present(root)) then
      mype_ = myid
      root_ = root
      if (myid /= root) verbose = .false.
   endif

   gsi_ = .false.
   if (present(gsiset)) gsi_ = gsiset

   ! Get dimensions from file
   call read_dims_(fname, nlat_, nlon_, nlev_, rc, mype_, root_)
   if (rc /= 0) return

   init_ = bvars%initialized
   if (init_) then
      nlat = bvars%nlat
      nlon = bvars%nlon
      nlev = bvars%nsig
      if (nlon_ /= nlon .or. nlat_ /= nlat .or. nlev_ /= nlev) then
         if (mype_ == root_) then
            print *, myname_, ': nlat(file)/nlat(req) = ', nlat_, nlat
            print *, myname_, ': nlon(file)/nlon(req) = ', nlon_, nlon
            print *, myname_, ': nlev(file)/nlev(req) = ', nlev_, nlev
            print *, myname_, ': inconsistent dimensions, aborting'
         endif
         rc = 1
         return
      endif
   else
      nlat = nlat_
      nlon = nlon_
      nlev = nlev_
      call init_GFSens_vars_(bvars, nlon, nlat, nlev, gsi=gsi_)
   endif

   call check_(nf90_open(fname, NF90_NOWRITE, ncid), rc, mype_, root_)
   if (rc /= 0) return

   ! Read 3D variables
   allocate(data3d(nlon, nlat, nlev))
   do nv = 1, bvars%nv3d
      fvname = gfs_varname_(bvars%gsi_vnames3d(nv))
      ier = nf90_inq_varid(ncid, trim(fvname), varid)
      if (ier /= NF90_NOERR) then
         if (verbose) print *, myname_, ': variable not found in file: ', trim(fvname), &
              ' (GSI name: ', trim(bvars%gsi_vnames3d(nv)), ')'
         data3d = 0.0
      else
         call check_(nf90_get_var(ncid, varid, data3d), rc, mype_, root_)
      endif
      if (gsi_) then
         do kk = 1, nlev
            bvars%ptr3d(:,:,kk,nv) = transpose(data3d(:,:,kk))
         enddo
      else
         bvars%ptr3d(:,:,:,nv) = data3d
      endif
   enddo
   deallocate(data3d)

   ! Read 2D variables
   allocate(data2d(nlon, nlat))
   do nv = 1, bvars%nv2d
      fvname = gfs_varname_(bvars%gsi_vnames2d(nv))
      ier = nf90_inq_varid(ncid, trim(fvname), varid)
      if (ier /= NF90_NOERR) then
         if (verbose) print *, myname_, ': variable not found in file: ', trim(fvname), &
              ' (GSI name: ', trim(bvars%gsi_vnames2d(nv)), ')'
         data2d = 0.0
      else
         call check_(nf90_get_var(ncid, varid, data2d), rc, mype_, root_)
      endif
      if (gsi_) then
         bvars%ptr2d(:,:,nv) = transpose(data2d)
      else
         bvars%ptr2d(:,:,nv) = data2d
      endif
   enddo
   deallocate(data2d)

   call check_(nf90_close(ncid), rc, mype_, root_)

   if (verbose) print *, myname_, ': finished reading file: ', trim(fname)

   ! Convert from GFS file units to GSI units and flip orientation if needed
   call gfs2gsi_(bvars)

end subroutine read_GFSens_

!---------------------------------------------------------------------------
! Convert GFS file units to GSI units.
! GFS NetCDF4 files (CF-compliant) are assumed to have:
!   - latitude south-to-north (same as GSI convention)
!   - longitude 0-360 eastward (same as GSI convention)
!   - vertical levels from model top to near-surface (same as GSI convention)
! Therefore only unit conversions are applied, not geometric transformations.
!
! Unit conversions:
!   - Surface pressure: Pa -> centibars (1 cb = 1 kPa = 1000 Pa)
subroutine gfs2gsi_(x)
   implicit none
   type(nc_GFSens_vars), intent(inout) :: x
   integer :: id

   ! Surface pressure: Pa -> centibars (1 cb = 1000 Pa)
   id = getindex(x%gsi_vnames2d, 'ps')
   if (id > 0) x%ptr2d(:,:,id) = x%ptr2d(:,:,id) * Pa_to_cb

end subroutine gfs2gsi_

!---------------------------------------------------------------------------
! flip_ and the associated latflip/levflip subroutines are provided for cases
! where GFS files do NOT follow the CF-compliant ordering assumed by gfs2gsi_.
! If a file has latitude going north-to-south or vertical levels ordered
! bottom-to-top, callers may invoke flip_() after calling gfs2gsi_().
! In the standard GFS NetCDF4 (FV3) gaussian grid output, flip_ is not needed.
subroutine flip_(x)
   implicit none
   type(nc_GFSens_vars), intent(inout) :: x
   integer :: im, jm, km, nv

   im = x%nlon
   jm = x%nlat
   km = x%nsig

   do nv = 1, x%nv2d
      call latflip2_(x%ptr2d(:,:,nv), im, jm, x%gsiset)
   enddo
   do nv = 1, x%nv3d
      call latflip3_(x%ptr3d(:,:,:,nv), im, jm, km, x%gsiset)
      call levflip_ (x%ptr3d(:,:,:,nv), im, jm, km, x%gsiset)
   enddo
end subroutine flip_

!---------------------------------------------------------------------------
! Flip latitude dimension (N->S to S->N).
! In GSI orientation (gsiset=.true.): array is (nlat,nlon[,nlev])
! In file orientation (gsiset=.false.): array is (nlon,nlat[,nlev])
subroutine latflip2_(q, im, jm, gsi)
   implicit none
   integer, intent(in) :: im, jm
   logical, intent(in) :: gsi
   real(4), intent(inout) :: q(:,:)
   real(4), allocatable :: dum(:)
   integer :: i, j

   if (gsi) then
      ! Array is (nlat, nlon): flip first dimension
      allocate(dum(im))
      do i = 1, jm/2
         dum(:) = q(i,:)
         q(i,:) = q(jm+1-i,:)
         q(jm+1-i,:) = dum(:)
      enddo
      deallocate(dum)
   else
      ! Array is (nlon, nlat): flip second dimension
      allocate(dum(im))
      do j = 1, jm/2
         dum(:) = q(:,j)
         q(:,j) = q(:,jm+1-j)
         q(:,jm+1-j) = dum(:)
      enddo
      deallocate(dum)
   endif
end subroutine latflip2_

!---------------------------------------------------------------------------
subroutine latflip3_(q, im, jm, km, gsi)
   implicit none
   integer, intent(in) :: im, jm, km
   logical, intent(in) :: gsi
   real(4), intent(inout) :: q(:,:,:)
   real(4), allocatable :: dum(:)
   integer :: i, j, k

   if (gsi) then
      ! Array is (nlat, nlon, nlev): flip first dimension
      allocate(dum(im))
      do k = 1, km
         do i = 1, jm/2
            dum(:) = q(i,:,k)
            q(i,:,k) = q(jm+1-i,:,k)
            q(jm+1-i,:,k) = dum(:)
         enddo
      enddo
      deallocate(dum)
   else
      ! Array is (nlon, nlat, nlev): flip second dimension
      allocate(dum(im))
      do k = 1, km
         do j = 1, jm/2
            dum(:) = q(:,j,k)
            q(:,j,k) = q(:,jm+1-j,k)
            q(:,jm+1-j,k) = dum(:)
         enddo
      enddo
      deallocate(dum)
   endif
end subroutine latflip3_

!---------------------------------------------------------------------------
! Flip vertical levels (top->bottom to bottom->top, i.e., surface at index 1 to nsig).
subroutine levflip_(q, im, jm, km, gsi)
   implicit none
   integer, intent(in) :: im, jm, km
   logical, intent(in) :: gsi
   real(4), intent(inout) :: q(:,:,:)
   real(4), allocatable :: dum(:,:)
   integer :: k

   if (gsi) then
      ! Array is (nlat, nlon, nlev)
      allocate(dum(jm, im))
      do k = 1, km/2
         dum(:,:) = q(:,:,k)
         q(:,:,k) = q(:,:,km+1-k)
         q(:,:,km+1-k) = dum(:,:)
      enddo
      deallocate(dum)
   else
      ! Array is (nlon, nlat, nlev)
      allocate(dum(im, jm))
      do k = 1, km/2
         dum(:,:) = q(:,:,k)
         q(:,:,k) = q(:,:,km+1-k)
         q(:,:,km+1-k) = dum(:,:)
      enddo
      deallocate(dum)
   endif
end subroutine levflip_

!---------------------------------------------------------------------------
subroutine get_pointer_2d_(vname, bvars, ptr, rc)
   implicit none
   character(len=*), intent(in) :: vname
   type(nc_GFSens_vars), intent(in) :: bvars
   real(4), pointer, intent(inout) :: ptr(:,:)
   integer, intent(out) :: rc
   integer :: id
   rc = -1
   id = getindex(bvars%gsi_vnames2d, trim(vname))
   if (id > 0) then
      ptr => bvars%ptr2d(:,:,id)
      rc = 0
   endif
end subroutine get_pointer_2d_

!---------------------------------------------------------------------------
subroutine get_pointer_3d_(vname, bvars, ptr, rc)
   implicit none
   character(len=*), intent(in) :: vname
   type(nc_GFSens_vars), intent(in) :: bvars
   real(4), pointer, intent(inout) :: ptr(:,:,:)
   integer, intent(out) :: rc
   integer :: id
   rc = -1
   id = getindex(bvars%gsi_vnames3d, trim(vname))
   if (id > 0) then
      ptr => bvars%ptr3d(:,:,:,id)
      rc = 0
   endif
end subroutine get_pointer_3d_

!---------------------------------------------------------------------------
subroutine check_(status, rc, myid, root)
   integer, intent(in)  :: status
   integer, intent(out) :: rc
   integer, intent(in)  :: myid, root
   rc = 0
   if (status /= nf90_noerr) then
      if (myid == root) print *, trim(nf90_strerror(status))
      rc = 999
   endif
end subroutine check_

end module m_nc_GFSens
