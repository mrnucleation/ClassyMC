# NaCl with P3M Electrostatics Examples

## NaCl_P3M - NVT Ensemble
Basic P3M example with fixed volume.

**Run:**
```bash
cd Example/NaCl_P3M
../../classyMC in.nacl
```

**Files:**
- `in.nacl` - Main input script
- `NaCl_P3M.clFF` - Forcefield definition with P3M parameters
- `NaCl_Config.clssy` - Initial configuration

**P3M Parameters in .clFF file:**
- `rcut 4.5` - Real-space cutoff (Angstroms)
- `alpha 0.35` - Screening parameter (1/Angstrom)
- `mesh 32 32 32` - FFT mesh dimensions
- `order 4` - Charge assignment order
- `precision 1.0E-5` - Target precision

## NaCl_P3M_NPT - NPT Ensemble
P3M with volume moves for constant pressure.

**Run:**
```bash
cd Example/NaCl_P3M_NPT
../../classyMC in.nacl_npt
```

**Additional Features:**
- Volume moves enabled (IsoVol)
- Pressure control
- Demonstrates P3M volume move capability

## Modifying Parameters

Edit the `.clFF` file to change P3M parameters:

```
forcefield 1
  rcut 4.5          # Change real-space cutoff
  alpha 0.35        # Adjust screening (≈ 3.5/rcut)
  mesh 32 32 32     # Change mesh size (powers of 2)
  order 4           # 2, 4, or 6
  precision 1.0E-5  # Tighter = more accurate
  1  1.0  0.5       # type charge rmin
  2 -1.0  0.5
end_forcefield
```

## Comparison with Ewald

Compare with `Example/NaCl_Ewald/` to validate P3M:
- Energies should agree to <0.1%
- P3M should be faster for this system size
- Both give same physical results

