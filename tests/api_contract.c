/* Small library-entry regression for the PNG transfer-function contract. */
#include <math.h>
#include <stdio.h>
#include "EZmock.h"
#include "errmsg.h"
#include "structs.h"

#define CHECK(test, message) do { \
  if (!(test)) { fprintf(stderr, "FAIL: %s\n", message); return 1; } \
} while (0)

int main(void) {
  int err = EZMOCK_SUCCESS;
  EZMOCK *odd = EZmock_init(1000, 15, 1, 1995, 1, &err);
  CHECK(!odd && err == EZMOCK_ERR_ARG_NGRID, "odd grid must be rejected");

  err = EZMOCK_SUCCESS;
  EZMOCK *ez = EZmock_init(1000, 16, 1, 1995, 1, &err);
  CHECK(ez && err == EZMOCK_SUCCESS, "initialize library");
  CHECK(EZmock_set_cosmology(ez, 1, 1, true, 1, 0, .3175, 0, -1, &err)
      == EZMOCK_SUCCESS, "set z=1 cosmology");
  EZMOCK_COSMO *cosmo = (EZMOCK_COSMO *) ez->cosmo;
  CHECK(cosmo->growth2 > .36 && cosmo->growth2 < .38,
      "growth2 must differ from one for this regression");

  const double k[] = {1e-5, .001, .01, .1, 1};
  const double transfer[] = {1, 1, 1, 1, 1};
  CHECK(EZmock_setup_linear_pk(ez, k, 5, transfer, NULL, 0, false, &err)
      == EZMOCK_SUCCESS, "set transfer function");
  EZMOCK_PK *pk = (EZMOCK_PK *) ez->pk;
  for (int i = 0; i < pk->n; i++)
    CHECK(pk->P[i] == 1, "library must not multiply T(k) by growth2");
  CHECK(EZmock_setup_linear_pk(ez, k, 5, transfer, transfer, .1, false, &err)
      == EZMOCK_SUCCESS, "set transfer function with BAO option");
  pk = (EZMOCK_PK *) ez->pk;
  for (int i = 0; i < pk->n; i++)
    CHECK(pk->P[i] == 1, "BAO branch must not multiply T(k) by growth2");

  err = EZMOCK_SUCCESS;
  CHECK(EZmock_setup_linear_pk(ez, k, 5, transfer, NULL, 0, true, &err)
      == EZMOCK_ERR_PNG_UNSUPPORTED, "reject log interpolation");
  err = EZMOCK_SUCCESS;
  real noise = 0;
  CHECK(EZmock_create_dens_field(ez, NULL, false, &noise, false, false, &err)
      == EZMOCK_ERR_PNG_UNSUPPORTED, "reject legacy white-noise path");

  err = EZMOCK_SUCCESS;
  EZmock_set_fnl(ez, 100);
  EZmock_set_b_phi(ez, 2.65);
  CHECK(EZmock_create_dens_field(ez, NULL, false, NULL, true, false, &err)
      == EZMOCK_SUCCESS, "create tracer-injection density field");
  EZMOCK_MESH *mesh = (EZMOCK_MESH *) ez->mesh;
  const size_t ncell = (size_t) 16 * 16 * 16;
  double max_diff = 0;
  for (size_t i = 0; i < ncell; i++) {
    double diff = fabs(mesh->phi_png[i] - mesh->phi[i]);
    if (diff > max_diff) max_diff = diff;
  }
  CHECK(max_diff < 1e-12, "B_PHI must disable implicit field-level PNG");
  EZmock_destroy(ez);
  puts("PASS: API transfer, unsupported paths, and PNG mode selection");
  return 0;
}
