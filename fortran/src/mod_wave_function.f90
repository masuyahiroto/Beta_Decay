module mod_wave_function
  use mod_constants, only: MAX_SP, MAX_OBDM
  implicit none

  ! Single-particle state quantum numbers
  type :: SingleParticleState
    integer :: n    ! radial quantum number
    integer :: l    ! orbital angular momentum
    integer :: j2   ! 2 * total angular momentum
    integer :: tz2  ! 2 * isospin projection: -1 proton, +1 neutron
  end type SingleParticleState

  ! One-body density matrix entry
  type :: OBDMEntry
    integer  :: p_idx   ! index into sp_states (proton)
    integer  :: n_idx   ! index into sp_states (neutron)
    real(8)  :: val
  end type OBDMEntry

  ! Nuclear state
  type :: NuclearState
    integer :: J2, parity, Z, N
    integer :: n_sp                        ! number of sp states stored
    integer :: n_obdm                      ! number of OBDM entries
    type(SingleParticleState) :: sp(MAX_SP)
    type(OBDMEntry)           :: obdm(MAX_OBDM)
  end type NuclearState

contains

  ! Add a single-particle state, return its index
  subroutine add_sp_state(state, n, l, j2, tz2, idx)
    type(NuclearState), intent(inout) :: state
    integer, intent(in)  :: n, l, j2, tz2
    integer, intent(out) :: idx
    state%n_sp = state%n_sp + 1
    idx = state%n_sp
    state%sp(idx)%n   = n
    state%sp(idx)%l   = l
    state%sp(idx)%j2  = j2
    state%sp(idx)%tz2 = tz2
  end subroutine add_sp_state

  ! Set OBDM element by state indices
  subroutine set_obdm(state, p_idx, n_idx, val)
    type(NuclearState), intent(inout) :: state
    integer, intent(in) :: p_idx, n_idx
    real(8), intent(in) :: val
    integer :: k
    ! update if already exists
    do k = 1, state%n_obdm
      if (state%obdm(k)%p_idx == p_idx .and. &
          state%obdm(k)%n_idx == n_idx) then
        state%obdm(k)%val = val
        return
      end if
    end do
    state%n_obdm = state%n_obdm + 1
    state%obdm(state%n_obdm)%p_idx = p_idx
    state%obdm(state%n_obdm)%n_idx = n_idx
    state%obdm(state%n_obdm)%val   = val
  end subroutine set_obdm

  ! Get OBDM element (0 if not set)
  function get_obdm(state, p_idx, n_idx) result(val)
    type(NuclearState), intent(in) :: state
    integer, intent(in) :: p_idx, n_idx
    real(8) :: val
    integer :: k
    val = 0.0d0
    do k = 1, state%n_obdm
      if (state%obdm(k)%p_idx == p_idx .and. &
          state%obdm(k)%n_idx == n_idx) then
        val = state%obdm(k)%val
        return
      end if
    end do
  end function get_obdm

  ! Initialize a NuclearState to zero
  subroutine init_state(state, J2, parity, Z, N)
    type(NuclearState), intent(out) :: state
    integer, intent(in) :: J2, parity, Z, N
    state%J2     = J2
    state%parity = parity
    state%Z      = Z
    state%N      = N
    state%n_sp   = 0
    state%n_obdm = 0
  end subroutine init_state

end module mod_wave_function
