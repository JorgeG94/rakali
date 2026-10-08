!! Unstructured triangular mesh type for the central-upwind solver
module rki_mesh_unstr
   !! Cell-centric unstructured mesh, flat 1D layout.
   !!
   !! Stores per-cell geometry (centroid, area, CFL radius) and
   !! per-edge-of-cell quantities (length, outward normal, midpoint, neighbour
   !! cell, neighbour-side edge index, boundary group). Edges are implicit:
   !! every triangle has exactly three edges, indexed locally as i in 1..3.
   !! Edge i is the edge opposite local vertex i.
   !!
   !! Per-edge-of-cell arrays use shape (3, ncells_total) so the three
   !! edges of cell k are contiguous in memory and the cell index is the
   !! slow-varying (column-major) axis -- the natural Fortran layout for
   !! a do concurrent (k = 1:ncells_owned) hot loop.
   !!
   !! All cells are owned (`ncells_owned == ncells_total`). Geometry is built
   !! once at init and never modified during the run.
   use rki_constants, only: wp, RAD2DEG
   use pic_types, only: int64
   use pic_logger, only: logger => global_logger
   use pic_io, only: to_char
   implicit none
   private

   public :: mesh_unstr_t
   public :: mesh_unstr_init_from_triangles
   public :: mesh_unstr_destroy
   public :: mesh_unstr_enter_data
   public :: mesh_unstr_exit_data
   public :: mesh_unstr_si_quality_gate

   integer, parameter :: EDGES_PER_TRI = 3
      !! Triangles only.

   real(wp), parameter :: SI_MAX_DEV_ANGLE_DEG = 5.0_wp
      !! Fail-loud ceiling on the max circumcentre-line deviation angle for
      !! the semi-implicit TPFA operator.  A clean Delaunay mesh measures
      !! < 1e-2 deg (edge orthogonal to the circumcentre line by
      !! construction); anything above 5 deg means a non-Delaunay / corrupted
      !! triangulation that TPFA would solve inconsistently — abort at
      !! configure rather than silently return a subtly-wrong eta.

   type :: mesh_unstr_t
      !! Triangular mesh in flat layout.

      ! ----- sizes -----
      integer :: nverts = 0
         !! Number of unique vertices.
      integer :: ncells_owned = 0
         !! Number of cells (live at indices 1..ncells_owned).
      integer :: ncells_total = 0
         !! Equal to `ncells_owned`.
      integer :: n_boundary_edges = 0
         !! Total boundary edges across owned cells.

      ! ----- vertex arrays -----
      real(wp), allocatable :: vert_x(:)
         !! Vertex x coordinate, size nverts.
      real(wp), allocatable :: vert_y(:)
         !! Vertex y coordinate, size nverts.

      ! ----- per-cell arrays -----
      integer, allocatable :: tri_v(:, :)
         !! Vertex ids of each triangle, shape (3, ncells_total). Ordered
         !! counter-clockwise; local vertex i is opposite local edge i.
      real(wp), allocatable :: centroid_x(:)
         !! Cell centroid x, size ncells_total.
      real(wp), allocatable :: centroid_y(:)
         !! Cell centroid y, size ncells_total.
      real(wp), allocatable :: area(:)
         !! Cell area, size ncells_total.
      real(wp), allocatable :: radius(:)
         !! Inscribed-circle radius, size ncells_total. Used as the CFL
         !! length scale for CFL computation.
      integer, allocatable :: n_boundaries(:)
         !! Number of boundary edges on each cell, 0..3, size ncells_total.

      ! ----- per-edge-of-cell arrays, shape (3, ncells_total) -----
      real(wp), allocatable :: edge_len(:, :)
         !! Length of edge i of cell k.
      real(wp), allocatable :: edge_n1(:, :)
         !! Outward unit normal x-component for edge i of cell k.
      real(wp), allocatable :: edge_n2(:, :)
         !! Outward unit normal y-component for edge i of cell k.
      real(wp), allocatable :: edge_mid_x(:, :)
         !! Edge midpoint x.
      real(wp), allocatable :: edge_mid_y(:, :)
         !! Edge midpoint y.
      integer, allocatable :: neighbours(:, :)
         !! Neighbouring cell index across edge i of cell k, or -1 if the
         !! edge is on a physical boundary.
      integer, allocatable :: neighbour_edges(:, :)
         !! Local edge index (1..3) of the neighbour that matches edge i
         !! of cell k, or -1 on a physical boundary.
      integer, allocatable :: surrogate_neighbours(:, :)
         !! Same as neighbours, but k itself if the edge is a physical
         !! boundary. Used by the gradient stencil so the reconstruction
         !! kernel never branches on boundaries.
      integer, allocatable :: bc_group(:, :)
         !! Boundary group id on edge i of cell k, 0 on interior edges.

      ! ----- semi-implicit (TPFA) geometry, computed once at init -----
      ! The semi-implicit free-surface solver (`time_scheme="semi_implicit"`)
      ! places eta at the triangle CIRCUMCENTRE (NOT the centroid the KNP
      ! path uses): a Delaunay edge is orthogonal to the circumcentre line
      ! by construction, which makes the two-point (TPFA) edge flux
      ! consistent (Casulli & Walters 2000).  Centroids stay untouched for
      ! the KNP reconstruction; these fields are the SI-specific additions.
      real(wp), allocatable :: circumcentre_x(:)
         !! Cell circumcentre x, size ncells_total.  eta lives here on the
         !! SI path.  Computed from the three vertices in `compute_geometry`.
      real(wp), allocatable :: circumcentre_y(:)
         !! Cell circumcentre y, size ncells_total.
      real(wp), allocatable :: edge_delta(:, :)
         !! Circumcentre-to-circumcentre distance across edge i of cell k,
         !! shape (3, ncells_total) (m).  Interior edge: |cc(nbr) - cc(k)|.
         !! Boundary edge: circumcentre-to-edge-midpoint distance (the
         !! Dirichlet ghost sits at the mirror across the edge).  The TPFA
         !! transmissibility is `Hf * edge_len / edge_delta`.  Built in a
         !! second pass (`compute_si_geometry`) once neighbours are known.
      real(wp), allocatable :: edge_dev_angle(:, :)
         !! Deviation angle (degrees) between the outward edge-normal and the
         !! circumcentre line across edge i of cell k, shape (3, ncells_total).
         !! 0 = orthogonal (ideal Delaunay).  The mesh-quality-gate metric;
         !! boundary edges carry 0 (not meaningful for the Dirichlet ghost).
      real(wp) :: max_dev_angle = 0.0_wp
         !! Max interior `edge_dev_angle` over the owned cells (degrees) — the
         !! fail-loud quality-gate scalar (abort if > ~5 deg on the SI path).
      real(wp) :: min_edge_delta = 0.0_wp
         !! Smallest SIGNED circumcentre-line projection on the outward
         !! edge-normal over the owned interior edges (m) — the quality-gate
         !! scalar.  Non-positive = a degenerate (coincident circumcentres,
         !! right-triangle split) or locally non-Delaunay (flipped) pair =>
         !! zero / wrong-signed TPFA transmissibility => abort.

      integer(int64), allocatable :: cell_global_id(:)
         !! Cell ids `1..ncells_total`, written to the output / restart
         !! files.
   end type mesh_unstr_t

contains

   subroutine mesh_unstr_init_from_triangles(mesh, vert_x, vert_y, tri_v, bc_tags)
      !! Build a `mesh_unstr_t` from raw vertex coordinates and triangle
      !! vertex lists. Computes geometry, neighbour topology and the
      !! surrogate-neighbour table.
      !!
      !! `tri_v(1:3, k)` are the vertex ids (1-based) of triangle k in CCW
      !! order.
      !!
      !! `bc_tags(1:3, k)` is an optional per-edge boundary group id per
      !! cell; only edges that are physical boundaries (`neighbours == -1`)
      !! keep their tag, interior edges get 0.
      type(mesh_unstr_t), intent(out) :: mesh
      real(wp), intent(in) :: vert_x(:)
      real(wp), intent(in) :: vert_y(:)
      integer, intent(in) :: tri_v(:, :)
      integer, intent(in), optional :: bc_tags(:, :)

      integer :: n_owned, n_total

      if (size(tri_v, 1) /= EDGES_PER_TRI) then
         call logger%error("mesh_unstr_init_from_triangles: tri_v must have " &
                           //"3 rows, got "//to_char(size(tri_v, 1)))
         return
      end if

      n_total = size(tri_v, 2)
      n_owned = n_total

      mesh%nverts = size(vert_x)
      mesh%ncells_owned = n_owned
      mesh%ncells_total = n_total

      allocate (mesh%vert_x(mesh%nverts))
      allocate (mesh%vert_y(mesh%nverts))
      mesh%vert_x = vert_x
      mesh%vert_y = vert_y

      ! `tri_v` is the FULL list (owned + halo), so a direct copy.
      allocate (mesh%tri_v(EDGES_PER_TRI, n_total))
      mesh%tri_v = tri_v

      allocate (mesh%centroid_x(n_total), source=0.0_wp)
      allocate (mesh%centroid_y(n_total), source=0.0_wp)
      allocate (mesh%area(n_total), source=0.0_wp)
      allocate (mesh%radius(n_total), source=0.0_wp)
      allocate (mesh%n_boundaries(n_total), source=0)

      allocate (mesh%edge_len(EDGES_PER_TRI, n_total), source=0.0_wp)
      allocate (mesh%edge_n1(EDGES_PER_TRI, n_total), source=0.0_wp)
      allocate (mesh%edge_n2(EDGES_PER_TRI, n_total), source=0.0_wp)
      allocate (mesh%edge_mid_x(EDGES_PER_TRI, n_total), source=0.0_wp)
      allocate (mesh%edge_mid_y(EDGES_PER_TRI, n_total), source=0.0_wp)

      ! Semi-implicit (TPFA) geometry — circumcentres per cell; edge deltas
      ! + deviation angles filled by `compute_si_geometry` after neighbours.
      allocate (mesh%circumcentre_x(n_total), source=0.0_wp)
      allocate (mesh%circumcentre_y(n_total), source=0.0_wp)
      allocate (mesh%edge_delta(EDGES_PER_TRI, n_total), source=0.0_wp)
      allocate (mesh%edge_dev_angle(EDGES_PER_TRI, n_total), source=0.0_wp)

      allocate (mesh%neighbours(EDGES_PER_TRI, n_total), source=-1)
      allocate (mesh%neighbour_edges(EDGES_PER_TRI, n_total), source=-1)
      allocate (mesh%surrogate_neighbours(EDGES_PER_TRI, n_total), source=-1)
      allocate (mesh%bc_group(EDGES_PER_TRI, n_total), source=0)

      call compute_geometry(mesh)
      call build_neighbours(mesh)
      call build_surrogate_neighbours(mesh)
      call count_boundary_edges(mesh)
      ! Semi-implicit TPFA geometry: circumcentre-line edge deltas +
      ! deviation angles (needs neighbour circumcentres, so post-neighbours).
      call compute_si_geometry(mesh)

      if (present(bc_tags)) then
         block
            integer :: kk, ii
            if (size(bc_tags, 1) /= EDGES_PER_TRI .or. &
                size(bc_tags, 2) /= n_owned) then
               call logger%error("mesh_unstr_init_from_triangles: bc_tags "// &
                                 "shape mismatch, expected (3, ncells_owned)")
            else
               do kk = 1, n_owned
                  do ii = 1, EDGES_PER_TRI
                     if (mesh%neighbours(ii, kk) < 0) then
                        mesh%bc_group(ii, kk) = bc_tags(ii, kk)
                     end if
                  end do
               end do
            end if
         end block
      end if

      allocate (mesh%cell_global_id(n_total))
      block
         integer :: kk
         do kk = 1, n_total
            mesh%cell_global_id(kk) = int(kk, kind=int64)
         end do
      end block

      call logger%info("mesh_unstr: built "//to_char(n_owned)// &
                       " cells, "//to_char(mesh%n_boundary_edges)//" boundary edges")
   end subroutine mesh_unstr_init_from_triangles

   subroutine compute_geometry(mesh)
      !! Compute per-cell area, centroid, CFL radius, and per-edge
      !! length, outward normal, and midpoint for every cell -- owned
      !! AND halo. Halo cells need real geometry too because the
      !! reconstruction kernel reads neighbour centroids when computing
      !! edge values across partition edges. Edge i is opposite local
      !! vertex i; in CCW orientation that gives an outward normal.
      type(mesh_unstr_t), intent(inout) :: mesh

      integer :: k, i, va, vb, vc
      real(wp) :: x1, y1, x2, y2, x3, y3
      real(wp) :: ex, ey, len, perimeter, dmid, r_min

      do k = 1, mesh%ncells_total
         va = mesh%tri_v(1, k)
         vb = mesh%tri_v(2, k)
         vc = mesh%tri_v(3, k)
         x1 = mesh%vert_x(va)
         y1 = mesh%vert_y(va)
         x2 = mesh%vert_x(vb)
         y2 = mesh%vert_y(vb)
         x3 = mesh%vert_x(vc)
         y3 = mesh%vert_y(vc)

         ! Signed area via cross product; positive iff CCW.
         mesh%area(k) = 0.5_wp*((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1))
         mesh%centroid_x(k) = (x1 + x2 + x3)/3.0_wp
         mesh%centroid_y(k) = (y1 + y2 + y3)/3.0_wp

         ! Circumcentre (equidistant from the three vertices) — the SI
         ! free-surface places eta here so the TPFA edge flux is consistent
         ! on a Delaunay mesh (edge orthogonal to the circumcentre line).
         ! den = 2*(signed area)*2 = 4A; nonzero for a non-degenerate cell.
         block
            real(wp) :: den, a2, b2, c2
            den = 2.0_wp*(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2))
            if (abs(den) > tiny(1.0_wp)) then
               a2 = x1*x1 + y1*y1
               b2 = x2*x2 + y2*y2
               c2 = x3*x3 + y3*y3
               mesh%circumcentre_x(k) = &
                  (a2*(y2 - y3) + b2*(y3 - y1) + c2*(y1 - y2))/den
               mesh%circumcentre_y(k) = &
                  (a2*(x3 - x2) + b2*(x1 - x3) + c2*(x2 - x1))/den
            else
               ! Degenerate cell: fall back to the centroid (the quality
               ! gate will reject it via a non-positive edge_delta anyway).
               mesh%circumcentre_x(k) = mesh%centroid_x(k)
               mesh%circumcentre_y(k) = mesh%centroid_y(k)
            end if
         end block

         ! Edge i is opposite vertex i. Walking edges in CCW order:
         !   edge 1 = (v2 -> v3), edge 2 = (v3 -> v1), edge 3 = (v1 -> v2).
         ! Outward normal of an edge (xa,ya)->(xb,yb) in a CCW triangle is
         ! ( (yb - ya), -(xb - xa) ) / length.
         perimeter = 0.0_wp
         r_min = huge(1.0_wp)
         do i = 1, EDGES_PER_TRI
            block
               real(wp) :: xa, ya, xb, yb
               select case (i)
               case (1)
                  xa = x2
                  ya = y2
                  xb = x3
                  yb = y3
               case (2)
                  xa = x3
                  ya = y3
                  xb = x1
                  yb = y1
               case default
                  xa = x1
                  ya = y1
                  xb = x2
                  yb = y2
               end select
               ex = xb - xa
               ey = yb - ya
               len = sqrt(ex*ex + ey*ey)
               mesh%edge_len(i, k) = len
               mesh%edge_n1(i, k) = ey/len
               mesh%edge_n2(i, k) = -ex/len
               mesh%edge_mid_x(i, k) = 0.5_wp*(xa + xb)
               mesh%edge_mid_y(i, k) = 0.5_wp*(ya + yb)
               perimeter = perimeter + len
               ! Distance from the centroid to this edge's midpoint.
               dmid = sqrt((mesh%centroid_x(k) - mesh%edge_mid_x(i, k))**2 &
                           + (mesh%centroid_y(k) - mesh%edge_mid_y(i, k))**2)
               if (dmid < r_min) r_min = dmid
            end block
         end do

         ! CFL length scale: min centroid-to-edge-midpoint distance
         ! (= 1/3 of the shortest median). Equals the inscribed-circle
         ! radius for equilateral triangles, smaller (more conservative)
         ! for skewed ones.
         if (perimeter > 0.0_wp) then
            mesh%radius(k) = r_min
         end if
      end do
   end subroutine compute_geometry

   subroutine build_neighbours(mesh)
      !! Build the neighbours / neighbour_edges tables by matching edges
      !! that share the same unordered vertex pair. Uses a single sort of
      !! 64-bit composite keys instead of a hash table -- O(M log M) where
      !! M = 3 * ncells_total, and avoids dragging in any container code.
      !!
      !! Iterates over `ncells_total` so that edges between owned and
      !! halo cells (the partition edges) are linked just like interior
      !! edges. Halo↔halo matches are harmless: the kernel never walks
      !! into them because flux loops are bounded by `ncells_owned`.
      type(mesh_unstr_t), intent(inout) :: mesh

      integer(int64), allocatable :: keys(:)
      integer, allocatable :: cell_of(:), edge_of(:)
      integer(int64) :: nv64
      integer :: m, idx, k, i, va, vb, vmin, vmax
      integer :: p, q, k1, e1, k2, e2

      m = EDGES_PER_TRI*mesh%ncells_total
      if (m == 0) return

      allocate (keys(m), cell_of(m), edge_of(m))
      nv64 = int(mesh%nverts, kind=int64) + 1_int64

      idx = 0
      do k = 1, mesh%ncells_total
         do i = 1, EDGES_PER_TRI
            ! Edge i is opposite vertex i, so its endpoints are the other two.
            select case (i)
            case (1)
               va = mesh%tri_v(2, k)
               vb = mesh%tri_v(3, k)
            case (2)
               va = mesh%tri_v(3, k)
               vb = mesh%tri_v(1, k)
            case default
               va = mesh%tri_v(1, k)
               vb = mesh%tri_v(2, k)
            end select
            vmin = min(va, vb)
            vmax = max(va, vb)
            idx = idx + 1
            keys(idx) = int(vmin, kind=int64)*nv64 + int(vmax, kind=int64)
            cell_of(idx) = k
            edge_of(idx) = i
         end do
      end do

      call sort_by_key(keys, cell_of, edge_of)

      ! Walk consecutive runs of equal keys: a length-2 run is an interior
      ! edge (link both sides); length 1 is a physical boundary; length >2
      ! means a non-manifold mesh, which is a hard error here.
      p = 1
      do while (p <= m)
         q = p
         do while (q < m)
            if (keys(q + 1) /= keys(p)) exit
            q = q + 1
         end do
         select case (q - p + 1)
         case (1)
            ! boundary -- leave neighbours(:, k) at -1
         case (2)
            k1 = cell_of(p)
            e1 = edge_of(p)
            k2 = cell_of(p + 1)
            e2 = edge_of(p + 1)
            mesh%neighbours(e1, k1) = k2
            mesh%neighbour_edges(e1, k1) = e2
            mesh%neighbours(e2, k2) = k1
            mesh%neighbour_edges(e2, k2) = e1
         case default
            call logger%error("mesh_unstr: non-manifold edge shared by " &
                              //to_char(q - p + 1)//" triangles")
         end select
         p = q + 1
      end do

      deallocate (keys, cell_of, edge_of)
   end subroutine build_neighbours

   subroutine build_surrogate_neighbours(mesh)
      !! Surrogate neighbours: same as `neighbours` but the cell itself
      !! stands in for boundary edges. Lets the gradient
      !! kernel index neighbour values without branching on -1.
      !! Filled for every cell -- halo cells included -- so that any
      !! kernel that happens to dereference `surrogate_neighbours(:, k)`
      !! for a halo `k` still gets a valid index.
      type(mesh_unstr_t), intent(inout) :: mesh

      integer :: k, i

      do k = 1, mesh%ncells_total
         do i = 1, EDGES_PER_TRI
            if (mesh%neighbours(i, k) > 0) then
               mesh%surrogate_neighbours(i, k) = mesh%neighbours(i, k)
            else
               mesh%surrogate_neighbours(i, k) = k
            end if
         end do
      end do
   end subroutine build_surrogate_neighbours

   subroutine count_boundary_edges(mesh)
      !! Tally `n_boundaries(k)` for every cell (owned AND halo) and
      !! the total boundary-edge count across owned cells.
      !!
      !! The per-cell count is needed by the extrapolate kernel to
      !! branch into the right reconstruction path (0/1 boundary →
      !! three-neighbour, 2 → one-neighbour, 3 → constant). Halo
      !! cells need correct counts too because the extrapolate kernel
      !! runs over `ncells_total` so that ghost-cell edge values are
      !! available for the flux kernel to read at partition edges.
      type(mesh_unstr_t), intent(inout) :: mesh

      integer :: k, i, n

      mesh%n_boundary_edges = 0
      do k = 1, mesh%ncells_total
         n = 0
         do i = 1, EDGES_PER_TRI
            if (mesh%neighbours(i, k) < 0) n = n + 1
         end do
         mesh%n_boundaries(k) = n
         ! Only tally physical boundaries on owned cells for the
         ! mesh-wide count (halo boundaries are artefacts of the
         ! local sub-mesh, not real domain boundaries).
         if (k <= mesh%ncells_owned) then
            mesh%n_boundary_edges = mesh%n_boundary_edges + n
         end if
      end do
   end subroutine count_boundary_edges

   subroutine compute_si_geometry(mesh)
      !! Fill the semi-implicit TPFA edge geometry: the circumcentre-to-
      !! circumcentre distance `edge_delta` and the edge/centre-line
      !! deviation angle `edge_dev_angle` for every edge of every cell, plus
      !! the mesh-wide quality-gate scalars `max_dev_angle` / `min_edge_delta`
      !! over the OWNED interior edges.  Runs after `build_neighbours` so
      !! neighbour circumcentres are available.  Host-only (init-time; the
      !! device never re-derives geometry).
      !!
      !!   interior edge i of cell k, neighbour kn:
      !!     dvec = cc(kn) - cc(k);  edge_delta = |dvec|
      !!     dev  = acos( |n_hat . dvec| / |dvec| )   (0 = orthogonal, ideal)
      !!   boundary edge:
      !!     edge_delta = circumcentre-to-edge-midpoint distance (the
      !!     Dirichlet ghost sits at the mirror across the edge); dev = 0.
      !!
      !! The circumcentre line of a valid Delaunay edge is orthogonal to the
      !! edge by construction, so `dev` is ~0 everywhere on a clean mesh; a
      !! large `max_dev_angle` (or a non-positive `min_edge_delta` from an
      !! obtuse triangle whose circumcentre falls outside itself) flags a
      !! mesh the TPFA operator would solve inconsistently — the SI configure
      !! gate aborts on it (`mesh_unstr_si_quality_gate`).
      type(mesh_unstr_t), intent(inout) :: mesh

      integer :: k, i, kn
      real(wp) :: dxc, dyc, dlen, cosang, dmx, dmy
      real(wp) :: max_dev, min_delta

      max_dev = 0.0_wp
      min_delta = huge(1.0_wp)

      do k = 1, mesh%ncells_total
         do i = 1, EDGES_PER_TRI
            kn = mesh%neighbours(i, k)
            if (kn > 0) then
               dxc = mesh%circumcentre_x(kn) - mesh%circumcentre_x(k)
               dyc = mesh%circumcentre_y(kn) - mesh%circumcentre_y(k)
               dlen = sqrt(dxc*dxc + dyc*dyc)
               mesh%edge_delta(i, k) = dlen
               ! Deviation angle between the outward edge-normal and the
               ! centre line (both unit-normalised via the /dlen below).
               if (dlen > tiny(1.0_wp)) then
                  cosang = abs(mesh%edge_n1(i, k)*dxc + mesh%edge_n2(i, k)*dyc)/dlen
                  cosang = min(max(cosang, 0.0_wp), 1.0_wp)
                  mesh%edge_dev_angle(i, k) = acos(cosang)*RAD2DEG
               else
                  mesh%edge_dev_angle(i, k) = 0.0_wp
               end if
               ! Quality-gate reduction over OWNED interior edges only.
               ! min_delta tracks the SIGNED projection of the circumcentre
               ! line on the OUTWARD normal: both circumcentres sit on the
               ! shared edge's perpendicular bisector for ANY triangulation
               ! (so the deviation angle alone cannot flag a bad pair) — a
               ! locally NON-Delaunay pair puts the neighbour circumcentre
               ! on the WRONG side (negative projection => a flipped
               ! transmissibility sign), and a right-triangle pair puts
               ! them coincident (zero).  Both must abort.
               if (k <= mesh%ncells_owned) then
                  max_dev = max(max_dev, mesh%edge_dev_angle(i, k))
                  min_delta = min(min_delta, &
                                  dxc*mesh%edge_n1(i, k) + dyc*mesh%edge_n2(i, k))
               end if
            else
               ! Boundary (or partition-with-no-neighbour) edge: use the
               ! circumcentre-to-edge-midpoint distance; the Dirichlet ghost
               ! surface is imposed at the mirror across the edge.
               dmx = mesh%edge_mid_x(i, k) - mesh%circumcentre_x(k)
               dmy = mesh%edge_mid_y(i, k) - mesh%circumcentre_y(k)
               mesh%edge_delta(i, k) = sqrt(dmx*dmx + dmy*dmy)
               mesh%edge_dev_angle(i, k) = 0.0_wp
            end if
         end do
      end do

      if (min_delta == huge(1.0_wp)) min_delta = 0.0_wp
      mesh%max_dev_angle = max_dev
      mesh%min_edge_delta = min_delta
   end subroutine compute_si_geometry

   subroutine mesh_unstr_si_quality_gate(mesh, ok)
      !! Fail-loud mesh-quality gate for the semi-implicit TPFA operator
      !! (Casulli & Walters 2000 orthogonal-centre requirement, made
      !! explicit).  Aborts via `ok=.false.` (the caller error-stops) when
      !! either
      !!   * the max circumcentre-line deviation angle exceeds
      !!     `SI_MAX_DEV_ANGLE_DEG` (a non-Delaunay / corrupted mesh), or
      !!   * the smallest interior `edge_delta` is non-positive (an obtuse
      !!     triangle put its circumcentre outside itself => a flipped /
      !!     zero transmissibility `K_e = Hf*len/delta`).
      !! A clean Delaunay mesh (meshpy/pymetis pipeline) passes trivially.
      type(mesh_unstr_t), intent(in) :: mesh
      logical, intent(out) :: ok

      ok = .true.
      if (mesh%max_dev_angle > SI_MAX_DEV_ANGLE_DEG) then
         call logger%error("mesh_unstr: semi-implicit quality gate FAILED — "// &
                           "max circumcentre-line deviation angle "// &
                           to_char(mesh%max_dev_angle)//" deg exceeds the "// &
                           to_char(SI_MAX_DEV_ANGLE_DEG)//" deg ceiling.  The "// &
                           "TPFA free-surface operator is consistent only on a "// &
                           "near-Delaunay mesh (edge orthogonal to the "// &
                           "circumcentre line); this triangulation is not.")
         ok = .false.
      end if
      if (mesh%min_edge_delta <= 0.0_wp) then
         call logger%error("mesh_unstr: semi-implicit quality gate FAILED — "// &
                           "non-positive interior circumcentre-line distance "// &
                           "(min edge_delta = "//to_char(mesh%min_edge_delta)// &
                           " m).  An obtuse triangle has its circumcentre "// &
                           "outside itself => the TPFA transmissibility flips "// &
                           "sign / vanishes.  Use an acute / Delaunay mesh.")
         ok = .false.
      end if
   end subroutine mesh_unstr_si_quality_gate

   subroutine sort_by_key(keys, payload1, payload2)
      !! In-place quicksort keeping two integer payload arrays in lockstep
      !! with the int64 key array. Mesh-init only.
      integer(int64), intent(inout) :: keys(:)
      integer, intent(inout) :: payload1(:)
      integer, intent(inout) :: payload2(:)

      call quicksort(1, size(keys))
   contains
      recursive subroutine quicksort(lo, hi)
         integer, intent(in) :: lo, hi
         integer :: i, j
         integer(int64) :: pivot, tmpk
         integer :: tmp1, tmp2

         if (lo >= hi) return
         pivot = keys((lo + hi)/2)
         i = lo
         j = hi
         do
            do while (keys(i) < pivot)
               i = i + 1
            end do
            do while (keys(j) > pivot)
               j = j - 1
            end do
            if (i <= j) then
               tmpk = keys(i)
               keys(i) = keys(j)
               keys(j) = tmpk
               tmp1 = payload1(i)
               payload1(i) = payload1(j)
               payload1(j) = tmp1
               tmp2 = payload2(i)
               payload2(i) = payload2(j)
               payload2(j) = tmp2
               i = i + 1
               j = j - 1
            end if
            if (i > j) exit
         end do
         if (lo < j) call quicksort(lo, j)
         if (i < hi) call quicksort(i, hi)
      end subroutine quicksort
   end subroutine sort_by_key

   subroutine mesh_unstr_enter_data(mesh)
      !! Map the read-only mesh geometry and topology arrays into device
      !! memory. Called once at solver init alongside
      !! `state_unstr_enter_data`. Vertex coords + raw triangle vertex
      !! list stay host-only (kernels never read them after init).
      type(mesh_unstr_t), intent(inout) :: mesh

      ! Parent struct descriptor first (MOM6-miniapp pattern), then
      ! components attach to it. EVERY allocatable component must appear
      ! in a per-component map below -- an allocated-but-unattached field
      ! leaves an invalid device pointer in the parent descriptor and the
      ! first `do concurrent` deref of `mesh%<anything>` trips
      ! CUDA_ERROR_ILLEGAL_ADDRESS. Partition metadata no kernel reads
      ! still must be attached to keep the parent descriptor well-formed.
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
      ! [acc->omp-cpu] dropped device directive (host no-op): enter data
   end subroutine mesh_unstr_enter_data

   subroutine mesh_unstr_exit_data(mesh)
      !! Release the device copies of the mesh arrays. Counterpart to
      !! `mesh_unstr_enter_data`. Mesh data is read-only on the device,
      !! so nothing needs to be copied back.
      type(mesh_unstr_t), intent(inout) :: mesh

      ! Components first, parent struct last (mirror of enter_data).
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
      ! [acc->omp-cpu] dropped device directive (host no-op): exit data
   end subroutine mesh_unstr_exit_data

   subroutine mesh_unstr_destroy(mesh)
      !! Release all allocatable arrays.
      type(mesh_unstr_t), intent(inout) :: mesh

      if (allocated(mesh%vert_x)) deallocate (mesh%vert_x)
      if (allocated(mesh%vert_y)) deallocate (mesh%vert_y)
      if (allocated(mesh%tri_v)) deallocate (mesh%tri_v)
      if (allocated(mesh%centroid_x)) deallocate (mesh%centroid_x)
      if (allocated(mesh%centroid_y)) deallocate (mesh%centroid_y)
      if (allocated(mesh%area)) deallocate (mesh%area)
      if (allocated(mesh%radius)) deallocate (mesh%radius)
      if (allocated(mesh%n_boundaries)) deallocate (mesh%n_boundaries)
      if (allocated(mesh%edge_len)) deallocate (mesh%edge_len)
      if (allocated(mesh%edge_n1)) deallocate (mesh%edge_n1)
      if (allocated(mesh%edge_n2)) deallocate (mesh%edge_n2)
      if (allocated(mesh%edge_mid_x)) deallocate (mesh%edge_mid_x)
      if (allocated(mesh%edge_mid_y)) deallocate (mesh%edge_mid_y)
      if (allocated(mesh%neighbours)) deallocate (mesh%neighbours)
      if (allocated(mesh%neighbour_edges)) deallocate (mesh%neighbour_edges)
      if (allocated(mesh%surrogate_neighbours)) deallocate (mesh%surrogate_neighbours)
      if (allocated(mesh%bc_group)) deallocate (mesh%bc_group)
      if (allocated(mesh%circumcentre_x)) deallocate (mesh%circumcentre_x)
      if (allocated(mesh%circumcentre_y)) deallocate (mesh%circumcentre_y)
      if (allocated(mesh%edge_delta)) deallocate (mesh%edge_delta)
      if (allocated(mesh%edge_dev_angle)) deallocate (mesh%edge_dev_angle)

      if (allocated(mesh%cell_global_id)) deallocate (mesh%cell_global_id)

      mesh%nverts = 0
      mesh%ncells_owned = 0
      mesh%ncells_total = 0
      mesh%n_boundary_edges = 0
   end subroutine mesh_unstr_destroy

end module rki_mesh_unstr
