module mod_shape_factor
  implicit none

  ! Shape factor type identifiers
  integer, parameter :: SHAPE_ALLOWED          = 1  ! allowed GT
  integer, parameter :: SHAPE_FERMI            = 2  ! allowed Fermi
  integer, parameter :: SHAPE_FF_NONUNIQUE     = 3  ! first-forbidden non-unique
  integer, parameter :: SHAPE_FF_UNIQUE        = 4  ! first-forbidden unique

  ! Shape factor container (replaces Python closures)
  type :: ShapeFactor
    integer  :: type_id = SHAPE_ALLOWED
    real(8)  :: B_GT  = 0.0d0  ! for SHAPE_ALLOWED
    real(8)  :: B_F   = 0.0d0  ! for SHAPE_FERMI
    real(8)  :: xi    = 0.0d0  ! for SHAPE_FF_NONUNIQUE
    real(8)  :: zeta  = 0.0d0  ! for SHAPE_FF_NONUNIQUE
    real(8)  :: M_ff  = 0.0d0  ! for SHAPE_FF_UNIQUE
  end type ShapeFactor

contains

  ! Evaluate shape factor C(W) at electron energy W (in units of m_e c^2)
  function eval_shape_factor(sf, W) result(C)
    type(ShapeFactor), intent(in) :: sf
    real(8), intent(in) :: W
    real(8) :: C

    select case (sf%type_id)
    case (SHAPE_ALLOWED)
      C = sf%B_GT
    case (SHAPE_FERMI)
      C = sf%B_F
    case (SHAPE_FF_NONUNIQUE)
      ! C(W) = (xi - zeta*W)^2 + zeta^2*(W^2-1)/3
      C = (sf%xi - sf%zeta * W)**2 + sf%zeta**2 * (W**2 - 1.0d0) / 3.0d0
    case (SHAPE_FF_UNIQUE)
      C = sf%M_ff**2
    case default
      C = 0.0d0
    end select
  end function eval_shape_factor

  ! Constructors
  function make_allowed_sf(B_GT) result(sf)
    real(8), intent(in) :: B_GT
    type(ShapeFactor) :: sf
    sf%type_id = SHAPE_ALLOWED
    sf%B_GT    = B_GT
  end function make_allowed_sf

  function make_fermi_sf(B_F) result(sf)
    real(8), intent(in) :: B_F
    type(ShapeFactor) :: sf
    sf%type_id = SHAPE_FERMI
    sf%B_F     = B_F
  end function make_fermi_sf

  function make_ff_sf(xi, zeta) result(sf)
    real(8), intent(in) :: xi, zeta
    type(ShapeFactor) :: sf
    sf%type_id = SHAPE_FF_NONUNIQUE
    sf%xi      = xi
    sf%zeta    = zeta
  end function make_ff_sf

end module mod_shape_factor
