# NaCl P3M Large System - Energy Comparison

## System Details

- **Box**: 30 x 30 x 30 Angstroms
- **Ions**: 27 Na+ and 27 Cl- (54 total atoms)
- **Structure**: CsCl-like arrangement (ions on offset cubic grids)
  - Na+ at positions (10n, 10m, 10p) where n,m,p ∈ {0,1,2}
  - Cl- at positions (5+10n, 5+10m, 5+10p) where n,m,p ∈ {0,1,2}

## Files

| File | Description |
|------|-------------|
| `NaCl_Config.clssy` | ClassyMC configuration file |
| `NaCl_P3M.clFF` | ClassyMC P3M forcefield file |
| `in.nacl` | ClassyMC input script |
| `lammps_nacl.data` | LAMMPS data file (same configuration) |
| `in.lammps` | LAMMPS input script for energy comparison |

## Running ClassyMC

```bash
cd /home/troy/ClassyMC_OPUS/ClassyMC/Example/NaCl_P3M_Large
../../classyMC in.nacl
```

Look for the initial energy output:
```
P3M Mesh Energy: XXXX
P3M Self-energy Correction: YYYY
P3M Total Energy: ZZZZ kb
```

## Running LAMMPS

```bash
cd /home/troy/ClassyMC_OPUS/ClassyMC/Example/NaCl_P3M_Large
lmp -in in.lammps
```

LAMMPS will output:
```
Total potential energy: XXXX kcal/mol
Coulomb energy (short): YYYY kcal/mol  
Coulomb energy (long/PPPM): ZZZZ kcal/mol
```

## Unit Conversion

ClassyMC uses **kb** (Boltzmann constant) units.
LAMMPS real units use **kcal/mol**.

To convert:
```
E(kb) = E(kcal/mol) × 503.228
```

Or equivalently:
```
E(kcal/mol) = E(kb) / 503.228
```

## P3M Parameters (Both Codes)

| Parameter | ClassyMC | LAMMPS |
|-----------|----------|--------|
| Real-space cutoff | 10.0 Å | 10.0 Å |
| Alpha (screening) | 0.35 /Å | (auto) |
| Mesh size | 32×32×32 | (auto) |
| Precision | 1e-5 | 1e-5 |

Note: LAMMPS PPPM automatically determines optimal alpha and mesh based on the precision target.

## Expected Comparison

Both codes should give the same total electrostatic energy (within precision tolerance) when properly converted to the same units. Small differences may arise from:

1. Different optimal alpha values
2. Different mesh sizes
3. Different charge assignment schemes
4. Numerical precision

A difference of < 0.1% indicates good agreement.

## Troubleshooting

### ClassyMC gives different energy than LAMMPS

1. Check unit conversion (kb vs kcal/mol)
2. Verify configurations are identical
3. Check that rcut is the same
4. Ensure precision/accuracy settings match

### LAMMPS errors

If LAMMPS fails with "Unknown pair style", ensure LAMMPS is compiled with the KSPACE package:
```bash
make yes-kspace
make mpi
```

