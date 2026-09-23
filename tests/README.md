# PNG API contract check

The small C regression checks the public library entry at `z=1`: a constant
transfer function must remain constant even though `growth2` is not one. It
also checks rejection of unsupported input paths and the implicit `FNL_FIELD=0`
selection when `B_PHI` is enabled. It does not validate a complete mock power
spectrum or covariance.

Build the double-precision library first, then compile and run the check with
the same compiler and FFTW installation. For example, on NERSC with at most
four threads allocated:

```sh
make -C libEZmock TARGET=libEZmock.a USE_OMP=F
cc -std=c99 -IlibEZmock/EZmock -IlibEZmock/prand/src/header \
  -I/global/homes/l/lzy/anaconda3/envs/kirisame/include \
  tests/api_contract.c libEZmock.a \
  -L/global/homes/l/lzy/anaconda3/envs/kirisame/lib -lfftw3 -lm \
  -o api_contract
./api_contract
```
