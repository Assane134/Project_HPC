#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NX 100001
#define G 1.0
#define CFL 0.9
#define XLEFT -10.0
#define XRIGHT 10.0
#define XM 0.0
#define HLEFT 5.0
#define HRIGHT 2.0
#define TEND 3.0

// Initialisation des conditions (discontinues)
// We didn't use that function
void init_conditions(int nc, double *x, double *h, double *u) {
  for (int i = 0; i < nc; ++i) {
    if (x[i] < XM)
      h[i] = HLEFT;
    else
      h[i] = HRIGHT;
    u[i] = 0.0;
  }
}

// Modified Rusanov flux function for MPI
void rusanov_flux(int count, int rank, int size, double *h, double *u,
                  double flux[][2]) {

  // Each process computes its own flux
  // count is the number of cells for this rank
  // rank is the current process rank
  // size is the total number of processes

  int off =
      (rank > 0
           ? 0
           : 1); // Offset for the first process because it has no left neighbor
  int end = (rank == size - 1
                 ? count
                 : count + 1 - off); // Last process has no right neighbor
  for (int i = 0; i < end; ++i) {
    int i_off = i + off; // Adjust index for the first process to not include
                         // the first ghost cell at index 0
    double hL = h[i_off], hR = h[i_off + 1];
    double uL = u[i_off], uR = u[i_off + 1];

    double WL[2] = {hL, hL * uL};
    double WR[2] = {hR, hR * uR};

    double FL_L[2] = {hL * uL, hL * uL * uL + 0.5 * G * hL * hL};
    double FL_R[2] = {hR * uR, hR * uR * uR + 0.5 * G * hR * hR};

    double sL = fabs(uL) + sqrt(G * hL);
    double sR = fabs(uR) + sqrt(G * hR);
    double smax = fmax(sL, sR);

    for (int j = 0; j < 2; ++j)
      flux[i][j] = 0.5 * (FL_L[j] + FL_R[j]) - 0.5 * smax * (WR[j] - WL[j]);
  }
}

int main(int argc, char *argv[]) {
  int nx = (argc > 1 ? atoi(argv[1]) : NX);

  printf("Number of grid points: %d\n", nx);
  int rank, size;

  // Initialisation de MPI
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double t_start = MPI_Wtime(); // Start timing the execution

  printf("Rank %d of %d\n", rank, size);

  int nc = nx - 1;

  int chunk = nc / size; // Number of cells per process
  int rem = nc % size;   // Remainder cells
  printf("Chunk size: %d, Remainder: %d\n", chunk, rem);

  int start = rank * chunk; // Start index for this process
  int count = (rank == size - 1)
                  ? chunk + rem
                  : chunk; // number of cells for this process, we give the last
                           // process the remainder cells

  double dx = (XRIGHT - XLEFT) / (nx - 1);
  double x[nx], xc[nx - 1];

  int L = 2 + count; // 2 ghost cells + count cells for this process,
                     // we add 2 ghost cells to the left and right of the local
                     // domain

  // Allocate arrays for the current process
  double h[L], u[L], hn[L], un[L];
  double W[2][L], Wn[2][L]; // flux[2][NX - 1];
  double flux[L][2];

  // Grids and centers
  for (int i = 0; i < nx; ++i)
    x[i] = XLEFT + i * dx;
  for (int i = 0; i < nc; ++i)
    xc[i] = 0.5 * (x[i] + x[i + 1]);

  // Initialisation
  // We only initialize the cells that this process owns
  for (int i = 0; i < count; ++i) {
    int gi = i + 1; // we add 1 to skip the ghost cell at index 0
    double cellpos = XLEFT + (start + i) * dx;
    h[gi] = (cellpos < XM ? HLEFT : HRIGHT);
    u[gi] = 0.0;
    W[0][gi] = h[gi];
    W[1][gi] = h[gi] * u[gi];
  }

  double t = 0.0;

  int left = (rank == 0 ? MPI_PROC_NULL : rank - 1);         // Left neighbor
  int right = (rank == size - 1 ? MPI_PROC_NULL : rank + 1); // Right neighbor

  while (t < TEND) {
    // Vitesse u et célérité
    double local_max_speed = 0.0; // Local maximum speed for this process
    double max_speed = 0.0;       // Global maximum speed across all processes
    for (int i = 1; i <= count; ++i) {
      h[i] = W[0][i];
      u[i] = W[1][i] / (h[i] > 1e-6 ? h[i] : 1e-6);
      double c = sqrt(G * h[i]);
      double speed = fabs(u[i]) + c;
      if (speed > local_max_speed)
        local_max_speed = speed;
    }

    MPI_Allreduce(&local_max_speed, &max_speed, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD); // Get the maximum speed across all processes

    // Pas de temps
    double dt = CFL * dx / max_speed;
    if (t + dt > TEND)
      dt = TEND - t;
    double nu = dt / dx;
    t += dt;

    // Timing communication
    if (left != MPI_PROC_NULL) {
      MPI_Sendrecv(&W[0][1], 1, MPI_DOUBLE, left, 1, &W[0][0], 1, MPI_DOUBLE,
                   left, 1, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE); // Send state value at index 1 (true
                                       // value) to left neighbor
      MPI_Sendrecv(&W[1][1], 1, MPI_DOUBLE, left, 0, &W[1][0], 1, MPI_DOUBLE,
                   left, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE); // Send state value at index 1 (true
                                       // value) to left neighbor
    }

    if (right != MPI_PROC_NULL) {
      MPI_Sendrecv(&W[0][count], 1, MPI_DOUBLE, right, 1, &W[0][count + 1], 1,
                   MPI_DOUBLE, right, 1, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE); // Send state value at index count (true
                                       // value) to right neighbor
      MPI_Sendrecv(&W[1][count], 1, MPI_DOUBLE, right, 0, &W[1][count + 1], 1,
                   MPI_DOUBLE, right, 0, MPI_COMM_WORLD,
                   MPI_STATUS_IGNORE); // Send state value at index count (true
                                       // value) to right neighbor
    }

    if (rank > 0) {
      h[0] = W[0][0];
      u[0] = W[1][0] / (h[0] > 1e-6 ? h[0] : 1e-6);
    } // update ghost cell at index 0 using received data from left neighbor

    if (rank < size - 1) {
      h[count + 1] = W[0][count + 1];
      u[count + 1] =
          W[1][count + 1] / (h[count + 1] > 1e-6 ? h[count + 1] : 1e-6);
    } // update ghost cell at index count + 1 using received data from right
      // neighbor

    // Calcul des flux
    rusanov_flux(count, rank, size, h, u, flux);

    int off =
        (rank > 0 ? 0
                  : 1); // Offset for the first process to skip the ghost cell
                        // at index 0 and to skip the first cell of the local
                        // domain because we will use Neumann BC at that cell
    int i1 =
        (rank == size - 1
             ? count - 1
             : count - off); // Avoid updating the
                             // last ghost cell and the cell before the ghost
                             // cell because we will use Neumann BC at that cell
    for (int i = 1; i <= i1; ++i) {
      int i_off = i + off;
      Wn[0][i_off] = W[0][i_off] - nu * (flux[i][0] - flux[i - 1][0]);
      Wn[1][i_off] = W[1][i_off] - nu * (flux[i][1] - flux[i - 1][1]);
    }

    // Neumann boundary conditions
    if (rank == 0) {
      Wn[0][1] = Wn[0][2];
      Wn[1][1] = Wn[1][2];
    }

    if (rank == size - 1) {
      Wn[0][count] = Wn[0][count - 1];
      Wn[1][count] = Wn[1][count - 1];
    }

    // Copy your updated cells back into W
    for (int i = 1; i <= count; ++i) {
      W[0][i] = Wn[0][i];
      W[1][i] = Wn[1][i];
    }
    // Affichage terminal (optionnel)
    printf("t = %.4f\n", t);
  }

  // Gather results from all processes
  int *counts = malloc(size * sizeof(int));
  int *displs = malloc(size * sizeof(int));
  int off = 0;
  for (int r = 0; r < size; ++r) {
    int c = (r == size - 1 ? chunk + rem : chunk);
    counts[r] = c;
    displs[r] = off;
    off += c;
  }

  double *all_h = NULL;
  double *all_hu = NULL;
  if (rank == 0) {
    all_h = malloc(nc * sizeof(double));
    all_hu = malloc(nc * sizeof(double));
  }

  // Gather the "h" component (W[0]) and the "hu" component (W[1])
  MPI_Gatherv(&W[0][1], count, MPI_DOUBLE, all_h, counts, displs, MPI_DOUBLE, 0,
              MPI_COMM_WORLD);

  MPI_Gatherv(&W[1][1], count, MPI_DOUBLE, all_hu, counts, displs, MPI_DOUBLE,
              0, MPI_COMM_WORLD);

  // On rank 0, write out x, h, and u = (hu/h)
  if (rank == 0) {
    FILE *f = fopen("output_mpi.dat", "w");
    for (int i = 0; i < nc; ++i) {
      double h_i = all_h[i];
      double hu_i = all_hu[i];
      fprintf(f, "%f\t%f\t%f\n", xc[i], h_i, hu_i / h_i);
    }
    fclose(f);
    free(all_h); // Free the allocated memory for all_h and all_hu
    free(all_hu);
  }

  double t_total = MPI_Wtime() - t_start;

  if (rank == 0) {
    printf(" Total runtime:      %8.3f s\n", t_total); // Print total runtime
  }

  free(counts);
  free(displs);

  MPI_Finalize();
  return 0;
}
