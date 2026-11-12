# CP2K total-potential export prerequisite

The VP solver needs a total local effective potential on the same Cartesian grid as the exported source wavefunction. The CP2K 2024.2 workflow used for this project adds the input section `DFT/PRINT/NATIVE_XC_CUBE`; its source patch is `NATIVE_XC_CUBE.patch`. The patch exports the native local XC component. It does not by itself write the complete total potential; the data-preparation step combines this component with the Hartree and local ionic/pseudopotential terms.

The patch changes two CP2K source files (`src/input_cp2k_print_dft.F` and `src/qs_scf_post_gpw.F`). The first registers the new input section; the second evaluates the final-density local XC potential and writes it as a Cube file for ordinary GPW/PBE calculations. The data-preparation step converts the Cube from Hartree to eV and combines it with the other local potential components to produce `total_potential_eV.npy`.

Use only with the matching CP2K revision (the patch was developed for CP2K 2024.2):

```bash
cp -a /path/to/cp2k-2024.2 /path/to/cp2k-2024.2-vpstm
cd /path/to/cp2k-2024.2-vpstm
patch -p1 < /path/to/NATIVE_XC_CUBE.patch
make -j24 ARCH=local VERSION=psmp
```

Then add the following to a normal GPW/PBE input:

```text
&DFT
  &PRINT
    &NATIVE_XC_CUBE ON
      STRIDE 1 1 1
    &END NATIVE_XC_CUBE
  &END PRINT
&END DFT
```

Check the executable and output log for `NATIVE_XC_CUBE`, `NATIVE_XC`, SCF convergence, and normal termination. UKS runs produce separate alpha/beta XC cubes. The exported XC cube is in Hartree; convert to eV before combining it with the other potential terms. Do not apply this patch to another CP2K release without resolving source-context conflicts and rerunning a small validation calculation.
