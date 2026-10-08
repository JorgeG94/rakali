!! Coastal-regime MOM6-style console statistics.
module rki_coastal_console_stats
   !! The coastal analogue of `rki_ocean_console_stats`: reduce the coastal
   !! state's conservation + stability scalars (total mass, kinetic energy,
   !! max advective CFL, and — when multilayer — volume-mean salinity /
   !! temperature) and emit them through the SHARED `rki_console_stats`
   !! formatter, so the coastal console block looks identical to the ocean's.
   !!
   !! Reductions run on device (only the scalars cross the bus) and follow
   !! the interior-cell + NaN/dry guard pattern of `compute_total_mass`.
   !! Mass + layer sums reuse the existing `rki_init` helpers; KE and CFL
   !! add their own reductions here.
   use rki_constants, only: wp
   use rki_grid, only: hgrid_t
   use rki_state, only: state_t
   use rki_console_stats, only: console_stats_t, console_stats_report, &
                                conservation_budget_t
   use rki_init, only: compute_total_mass, compute_layer_sums
   implicit none
   private

   public :: coastal_console_stats_report
   public :: coastal_total_ke, coastal_max_cfl
      !! Exposed for analytical unit tests (test_console_stats).

contains

   subroutine coastal_console_stats_report(this, s, dt, t, step)
      !! Reduce the coastal totals and print the shared MOM6-style block.
      !! `this` latches the t=0 reference
      !! on the first call (see `rki_console_stats`).
      type(console_stats_t), intent(inout) :: this
      type(state_t), intent(inout) :: s
      real(wp), intent(in) :: dt, t
      integer, intent(in) :: step

      real(wp) :: total_mass, total_ke, max_cfl, cell_area
      real(wp) :: sum_h, sum_hS, sum_hT, mean_S, mean_T, total_salt, total_heat
      logical :: has_st
      type(conservation_budget_t) :: bud

      cell_area = s%grid%dx*s%grid%dy
      has_st = s%use_multilayer

      total_mass = compute_total_mass(s)
      ! One fused device sweep returns both the KE sum and the max CFL.
      call coastal_reduce_stats(s, dt, total_ke, max_cfl)

      ! Mass conservation budget: cumulative volume out through the domain
      ! boundaries (accumulated per RK stage in the solver from the flux
      ! divergence) + cumulative volume added by surface sources (precip,
      ! evap, sponge h-relaxation).  Both terms must be tracked for the
      ! closed budget to be valid; if either is absent, fall back to raw
      ! drift (mass_active=.false.).
      bud%mass_out = s%barotropic%mass_out
      bud%mass_src = s%barotropic%mass_src
      bud%mass_active = s%barotropic%mass_out_tracked &
                        .and. s%barotropic%mass_src_tracked

      sum_h = 0.0_wp
      sum_hS = 0.0_wp
      sum_hT = 0.0_wp
      if (s%use_multilayer) call compute_layer_sums(s, sum_h, sum_hS, sum_hT)

      if (s%use_multilayer .and. sum_h > 0.0_wp) then
         mean_S = sum_hS/sum_h
         mean_T = sum_hT/sum_h
      else
         mean_S = 0.0_wp
         mean_T = 0.0_wp
      end if
      ! Column-integrated salt / heat for the drift block (area-weighted so
      ! the absolute figure is Σ hS·dx·dy; the Error ratio is scale-free).
      total_salt = sum_hS*cell_area
      total_heat = sum_hT*cell_area

      call console_stats_report(this, t, step, total_mass, total_ke, &
                                mean_S, mean_T, max_cfl, total_salt, total_heat, &
                                has_salt=has_st, has_temp=has_st, &
                                mean_age=0.0_wp, has_age=.false., budget=bud)
   end subroutine coastal_console_stats_report

   ! -------- device-side reductions (interior cells, NaN/dry guarded) -----
   !
   ! One fused sweep per rank (barotropic 2D / multilayer 3D) reduces BOTH
   ! the KE sum and the max advective CFL, so the console cadence touches the
   ! grid once instead of twice.  The per-cell physics lives in the pure
   ! `ke_cell` / `cfl_cell` formula helpers (`!$acc routine seq`, dry/NaN
   ! guard baked in — they return 0, the identity for both `+` and `max`, so
   ! a guarded cell drops out of either reduction), shared across both ranks.
   ! Device residency is assumed (`present`, matching the ocean twin): the
   ! production caller runs inside `solver_enter_data`/`solver_exit_data`, and
   ! the unit test maps the state the same way.

   subroutine coastal_reduce_stats(s, dt, ke, mcfl)
      !! Fused KE-sum + max-CFL reduction, dispatched on the layer mode.
      type(state_t), intent(inout) :: s
      real(wp), intent(in) :: dt
      real(wp), intent(out) :: ke, mcfl
      real(wp) :: area, idx, idy
      area = s%grid%dx*s%grid%dy
      idx = 1.0_wp/s%grid%dx
      idy = 1.0_wp/s%grid%dy
      if (s%use_multilayer) then
         call reduce_multilayer(s%multilayer%h_layer, s%multilayer%hu_layer, &
                                s%multilayer%hv_layer, s%grid, s%multilayer%nz_ml, &
                                area, dt, idx, idy, ke, mcfl)
      else
         call reduce_barotropic(s%barotropic%h, s%barotropic%hu, s%barotropic%hv, &
                                s%grid, area, dt, idx, idy, ke, mcfl)
      end if
   end subroutine coastal_reduce_stats

   function coastal_total_ke(s) result(ke)
      !! Σ 0.5·(hu²+hv²)/h·area — barotropic (2D) or summed over layers.  Thin
      !! wrapper over the fused sweep (max-CFL half discarded); kept public for
      !! the analytical KE unit test.
      type(state_t), intent(inout) :: s
      real(wp) :: ke
      real(wp) :: mcfl
      call coastal_reduce_stats(s, 1.0_wp, ke, mcfl)
   end function coastal_total_ke

   function coastal_max_cfl(s, dt) result(mcfl)
      !! max (|u|·dt/dx + |v|·dt/dy) over interior wet cells (advective CFL,
      !! matching the ocean's definition).  Thin wrapper over the fused sweep
      !! (KE half discarded); kept public for the analytical CFL unit test.
      type(state_t), intent(inout) :: s
      real(wp), intent(in) :: dt
      real(wp) :: mcfl
      real(wp) :: ke
      call coastal_reduce_stats(s, dt, ke, mcfl)
   end function coastal_max_cfl

   subroutine reduce_barotropic(h, hu, hv, grid, area, dt, idx, idy, ke, mcfl)
      !! Single fused sweep over the 2D interior.
      type(hgrid_t), intent(in) :: grid
      real(wp), intent(in) :: h(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hu(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: hv(grid%nx_total, grid%ny_total)
      real(wp), intent(in) :: area, dt, idx, idy
      real(wp), intent(out) :: ke, mcfl
      integer :: i, j, i_lo, i_hi, j_lo, j_hi
      real(wp) :: ke_acc, cfl_acc
      call interior_bounds(grid, i_lo, i_hi, j_lo, j_hi)
      ke_acc = 0.0_wp
      cfl_acc = 0.0_wp
      do concurrent(j=j_lo:j_hi, i=i_lo:i_hi) reduce(+:ke_acc) reduce(max:cfl_acc)
         ke_acc = ke_acc + ke_cell(h(i, j), hu(i, j), hv(i, j), area)
         cfl_acc = max(cfl_acc, cfl_cell(h(i, j), hu(i, j), hv(i, j), dt, idx, idy))
      end do
      ke = ke_acc
      mcfl = cfl_acc
   end subroutine reduce_barotropic

   subroutine reduce_multilayer(h_layer, hu_layer, hv_layer, grid, nz, &
                                area, dt, idx, idy, ke, mcfl)
      !! Single fused sweep over the 3D interior (summed over layers).
      type(hgrid_t), intent(in) :: grid
      integer, intent(in) :: nz
      real(wp), intent(in) :: h_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hu_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: hv_layer(grid%nx_total, grid%ny_total, nz)
      real(wp), intent(in) :: area, dt, idx, idy
      real(wp), intent(out) :: ke, mcfl
      integer :: i, j, k, i_lo, i_hi, j_lo, j_hi
      real(wp) :: ke_acc, cfl_acc
      call interior_bounds(grid, i_lo, i_hi, j_lo, j_hi)
      ke_acc = 0.0_wp
      cfl_acc = 0.0_wp
      do concurrent(k=1:nz, j=j_lo:j_hi, i=i_lo:i_hi) reduce(+:ke_acc) reduce(max:cfl_acc)
         ke_acc = ke_acc + ke_cell(h_layer(i, j, k), hu_layer(i, j, k), &
                                   hv_layer(i, j, k), area)
         cfl_acc = max(cfl_acc, cfl_cell(h_layer(i, j, k), hu_layer(i, j, k), &
                                         hv_layer(i, j, k), dt, idx, idy))
      end do
      ke = ke_acc
      mcfl = cfl_acc
   end subroutine reduce_multilayer

   pure subroutine interior_bounds(grid, i_lo, i_hi, j_lo, j_hi)
      !! Physical-cell index range (ghost rows excluded) shared by both sweeps.
      type(hgrid_t), intent(in) :: grid
      integer, intent(out) :: i_lo, i_hi, j_lo, j_hi
      i_lo = grid%nghost + 1
      i_hi = grid%nx_total - grid%nghost
      j_lo = grid%nghost + 1
      j_hi = grid%ny_total - grid%nghost
   end subroutine interior_bounds

   pure function ke_cell(h, hu, hv, area) result(v)
      !! Per-cell kinetic energy 0.5·(hu²+hv²)/h·area; 0 on a dry / NaN-guarded
      !! cell (the additive identity, so it drops out of the KE sum).
      !$acc routine seq
      real(wp), intent(in) :: h, hu, hv, area
      real(wp) :: v
      if (h > 0.0_wp .and. h < 1.0e6_wp) then
         v = 0.5_wp*(hu*hu + hv*hv)/h*area
      else
         v = 0.0_wp
      end if
   end function ke_cell

   pure function cfl_cell(h, hu, hv, dt, idx, idy) result(v)
      !! Per-cell advective CFL |u|·dt·idx + |v|·dt·idy; 0 on a dry / NaN-guarded
      !! cell (≤ any real CFL, so it drops out of the max).
      !$acc routine seq
      real(wp), intent(in) :: h, hu, hv, dt, idx, idy
      real(wp) :: v
      if (h > 0.0_wp .and. h < 1.0e6_wp) then
         v = abs(hu)/h*dt*idx + abs(hv)/h*dt*idy
      else
         v = 0.0_wp
      end if
   end function cfl_cell

end module rki_coastal_console_stats
