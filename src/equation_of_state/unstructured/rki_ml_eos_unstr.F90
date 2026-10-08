!! Equation of state and reference-density profile for the unstructured
!! multilayer solver.
module rki_ml_eos_unstr
   !! Linear two-tracer EOS:
   !!   rho_k = rho_0 + beta_S*(S_k - S_ref) - alpha_T*(T_k - T_ref)
   !! Also builds the horizontally-uniform reference-density profile
   !! `rho_ref(z)`, subtracted by both BPG kernels before differencing to
   !! cut the residual sigma pressure-gradient error in stratified water.
   use rki_constants, only: wp, DRY_TOLERANCE, H_VANISHED
   use rki_state_unstr, only: unstr_state_t
   use rki_mesh_unstr, only: mesh_unstr_t
   use rki_ml_common_unstr, only: ML_MIN_DEPTH_UNSTR
   use rki_eos, only: eos_t, eos_density_point
   implicit none
   private

   public :: ml_equation_of_state_unstr
   public :: ml_equation_of_state_shared_unstr
   public :: ml_compute_reference_density_unstr
   public :: rho_ref_at
#ifdef RAKALI_BENCHMARK_API
   ! Internal helper, exposed for benchmarks/bench_eos.F90 (no full unstr_state_t).
   public :: ml_equation_of_state_impl
#endif

contains

   subroutine ml_equation_of_state_unstr(state)
      !! Linear two-tracer EOS, per-cell, per-layer:
      !!   rho_k = rho_0 + beta_S*(S_k - S_ref) - alpha_T*(T_k - T_ref)
      !! with S_k = hS_k/h_k, T_k = hT_k/h_k. Dry layers fall back to the
      !! reference state. Computed over all cells (interior + halo) so the
      !! BPG gradient kernel reads neighbour rho without a separate halo.
      !!
      !! Outer shim: dereferences `tracers(idx)%hTr` on the host and passes
      !! flat 2D arrays into the impl — the array-of-derived-types descriptor
      !! can't be walked on the device under NVHPC stdpar.
      type(unstr_state_t), intent(inout) :: state

      if (.not. state%use_multilayer) return

      if (state%multilayer%eos_use_shared) then
         ! Opt-in shared nonlinear EOS (`&tracer_nml eos_form`).  All
         ! unstructured EOS consumers reach the EOS through this one
         ! entry, so the dispatch here covers the explicit RK2 path,
         ! the semi-implicit path, and solver init alike.
         call ml_equation_of_state_shared_unstr(state%multilayer%h_layer, &
                                                state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                                state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                                state%multilayer%rho_layer, &
                                                state%multilayer%eos, &
                                                state%multilayer%nz_ml, state%ncells_total)
      else
         call ml_equation_of_state_impl(state%multilayer%h_layer, &
                                        state%multilayer%tracers(state%multilayer%idx_salinity)%hTr, &
                                        state%multilayer%tracers(state%multilayer%idx_temperature)%hTr, &
                                        state%multilayer%rho_layer, &
                                        state%multilayer%S_ref, state%multilayer%T_ref, &
                                        state%multilayer%rho_0, state%multilayer%beta_S, state%multilayer%alpha_T, &
                                        state%multilayer%nz_ml, state%ncells_total)
      end if
   end subroutine ml_equation_of_state_unstr

   pure subroutine ml_equation_of_state_shared_unstr(h_layer, hS_layer, hT_layer, &
                                                     rho_layer, eos, nz, n)
      !! Opt-in SHARED (ocean) equation of state on the unstructured
      !! `(nz, ncells)` layout — the axis order is opposite to the
      !! structured shim, so this evaluates `eos_density_point` per
      !! (layer, cell) instead of forwarding to the structured 3D impls.
      !! `eos_density_point` is `pure` + `!$acc routine seq` and takes
      !! the flat-POD `eos_t` by value, so the handle rides into
      !! the kernel as a register copy — no device mapping.  Evaluation
      !! is at `eos%p_ref` (default 0 ⇒ surface/potential density).
      !!
      !! Vanishing-layer gate is the shared path's `H_VANISHED` (not the
      !! legacy kernel's `DRY_TOLERANCE`), matching the structured shared
      !! impls: layers at or below it fall back to `eos%rho0`.
      integer, intent(in) :: nz, n
      real(wp), intent(in) :: h_layer(nz, n), hS_layer(nz, n), hT_layer(nz, n)
      real(wp), intent(out) :: rho_layer(nz, n)
      type(eos_t), intent(in) :: eos

      integer :: cell, k
      real(wp) :: h_k, S_k, T_k, inv_h

      do concurrent(cell=1:n, k=1:nz) local(h_k, S_k, T_k, inv_h)
         h_k = h_layer(k, cell)
         if (h_k > H_VANISHED) then
            inv_h = 1.0_wp/h_k
            S_k = hS_layer(k, cell)*inv_h
            T_k = hT_layer(k, cell)*inv_h
            rho_layer(k, cell) = eos_density_point(eos, T_k, S_k, eos%p_ref)
         else
            rho_layer(k, cell) = eos%rho0
         end if
      end do
   end subroutine ml_equation_of_state_shared_unstr

   pure subroutine ml_equation_of_state_impl(h_layer, hS_layer, hT_layer, rho_layer, &
                                             S_ref, T_ref, rho_0, beta_S, alpha_T, &
                                             nz, n)
      !! Inner kernel — flat-array dummies so the device-side
      !! `do concurrent` doesn't chase the tracer-registry indirection.
      integer, intent(in)  :: nz, n
      real(wp), intent(in)  :: h_layer(nz, n), hS_layer(nz, n), hT_layer(nz, n)
      real(wp), intent(out) :: rho_layer(nz, n)
      real(wp), intent(in)  :: S_ref, T_ref, rho_0, beta_S, alpha_T

      integer :: cell, k
      real(wp) :: h_k, S_k, T_k, inv_h

      do concurrent(cell=1:n, k=1:nz) local(h_k, S_k, T_k, inv_h)
         h_k = h_layer(k, cell)
         if (h_k > DRY_TOLERANCE) then
            inv_h = 1.0_wp/h_k
            S_k = hS_layer(k, cell)*inv_h
            T_k = hT_layer(k, cell)*inv_h
         else
            S_k = S_ref
            T_k = T_ref
         end if
         rho_layer(k, cell) = rho_0 + beta_S*(S_k - S_ref) - alpha_T*(T_k - T_ref)
      end do
   end subroutine ml_equation_of_state_impl

   pure function rho_ref_at(z, rho_ref, dz, nz_ref, rho_0) result(rho)
      !$omp declare target
      !! Linear interpolation of the reference-density profile at depth `z`.
      !! Out-of-range values fall back to the nearest endpoint.
      !! When the profile is unset (`dz == 0`) returns `rho_0` so the
      !! kernel reduces to its pre-B2 behaviour.
      real(wp), intent(in) :: z, dz, rho_0
      integer, intent(in) :: nz_ref
      real(wp), intent(in) :: rho_ref(0:nz_ref)
      real(wp) :: rho, zr, frac
      integer :: i

      if (dz <= 0.0_wp .or. nz_ref < 1) then
         rho = rho_0
         return
      end if

      zr = max(0.0_wp, min(z, dz*real(nz_ref, wp)))
      i = int(zr/dz)
      if (i >= nz_ref) then
         rho = rho_ref(nz_ref)
      else
         frac = zr/dz - real(i, wp)
         rho = (1.0_wp - frac)*rho_ref(i) + frac*rho_ref(i + 1)
      end if
   end function rho_ref_at

   subroutine ml_compute_reference_density_unstr(state, mesh)
      !! Build the horizontally-uniform reference-density profile
      !! `rho_ref(z)` from the current `rho_layer`, sampled on a uniform
      !! z-grid `0..rho_ref_dz*nz_ref` and area-weighted across owned
      !! cells. Subtracted by both BPG kernels so the integrand is the
      !! anomaly `rho - rho_ref(z)`, scaling the residual sigma PGE down.
      !! Called once at solver init after the first EOS evaluation, on
      !! the host (one-time D2H+H2D).
      type(unstr_state_t), intent(inout) :: state
      type(mesh_unstr_t), intent(in) :: mesh

      integer :: cell, k, iz, nz, nz_ref
      real(wp) :: H_max, H_col, dz, z, z_top, z_bot, rho_at_z, area_w, sum_w
      real(wp), allocatable :: sum_rho_w(:)

      if (.not. state%use_multilayer) return
      if (.not. allocated(state%multilayer%rho_ref)) return

      nz = state%multilayer%nz_ml
      nz_ref = state%multilayer%nz_ref

      ! Pull layer data back so we can scan it on the host.
      !$omp target update from(state%multilayer%h_layer, state%multilayer%rho_layer)

      H_max = 0.0_wp
      do cell = 1, mesh%ncells_owned
         H_col = 0.0_wp
         do k = 1, nz
            H_col = H_col + state%multilayer%h_layer(k, cell)
         end do
         if (H_col > H_max) H_max = H_col
      end do

      if (H_max <= 0.0_wp) then
         state%multilayer%rho_ref = state%multilayer%rho_0
         state%multilayer%rho_ref_dz = 0.0_wp
         !$omp target update to(state%multilayer%rho_ref)
         return
      end if

      ! Small margin so all columns fit strictly inside the grid even
      ! after stage variations.
      H_max = 1.05_wp*H_max
      dz = H_max/real(nz_ref, wp)
      state%multilayer%rho_ref_dz = dz

      allocate (sum_rho_w(0:nz_ref), source=0.0_wp)
      sum_w = 0.0_wp

      do cell = 1, mesh%ncells_owned
         H_col = 0.0_wp
         do k = 1, nz
            H_col = H_col + state%multilayer%h_layer(k, cell)
         end do
         if (H_col < ML_MIN_DEPTH_UNSTR) cycle

         area_w = mesh%area(cell)
         sum_w = sum_w + area_w

         ! Walk the column from the surface (k = nz) downward to the
         ! bed (k = 1), tracking the depth `z_top` of the upper layer
         ! interface. For each ref grid point, find the layer that
         ! contains `z` and use its density.
         do iz = 0, nz_ref
            z = real(iz, wp)*dz
            z_top = 0.0_wp
            rho_at_z = state%multilayer%rho_layer(nz, cell)
            if (z >= H_col) then
               rho_at_z = state%multilayer%rho_layer(1, cell)
            else
               do k = nz, 1, -1
                  z_bot = z_top + state%multilayer%h_layer(k, cell)
                  if (z >= z_top .and. z <= z_bot) then
                     rho_at_z = state%multilayer%rho_layer(k, cell)
                     exit
                  end if
                  z_top = z_bot
               end do
            end if
            sum_rho_w(iz) = sum_rho_w(iz) + area_w*rho_at_z
         end do
      end do

      if (sum_w > 0.0_wp) then
         do iz = 0, nz_ref
            state%multilayer%rho_ref(iz) = sum_rho_w(iz)/sum_w
         end do
      else
         state%multilayer%rho_ref = state%multilayer%rho_0
      end if

      deallocate (sum_rho_w)

      !$omp target update to(state%multilayer%rho_ref)
   end subroutine ml_compute_reference_density_unstr

end module rki_ml_eos_unstr
