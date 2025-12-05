# Troubleshooting OpenFOAM Dimensional Errors in compInterFoam

## Common Error Messages

### Error 1: Wrong Dimensions
```
--> FOAM FATAL IO ERROR: (openfoam-2406)
Entry 'p_ref' in "controlDict.phaseChangeCoeffs" has dimensions [0 0 0 0 0 0 0]
but expected [1 -1 -2 0 0 0 0]
```

**Cause**: The parameter `p_ref` is defined without proper dimensional specification.

**Solution**: Add dimensional brackets before the value:
```cpp
// WRONG:
p_ref    101325;

// CORRECT:
p_ref    [1 -1 -2 0 0 0 0] 101325;  // Pa
```

### Error 2: Wrong Token Type
```
--> FOAM FATAL IO ERROR: (openfoam-2406)
Wrong token type - expected scalar value, found on line 64: punctuation '['
```

**Cause**: Incorrect syntax - likely missing the parameter name or having extra brackets.

**Solution**: Ensure proper format:
```cpp
// WRONG:
gasConstant [0 2 -2 -1 0 0 0] = 174;
gasConstant: [0 2 -2 -1 0 0 0] 174;

// CORRECT:
gasConstant [0 2 -2 -1 0 0 0] 174;
```

## Correct Format for phaseChangeCoeffs

```cpp
phaseChangeCoeffs
{
    model               clausius_clapeyron;

    // All physical quantities with dimensions must follow:
    // parameterName [kg m s K mol A cd] value;

    hf                  [0 2 -2 0 0 0 0] 9.1e6;    // J/kg
    gasConstant         [0 2 -2 -1 0 0 0] 174;     // J/(kg·K)
    Tsol                [0 0 0 1 0 0 0] 1941;      // K
    Tvap                [0 0 0 1 0 0 0] 3560;      // K
    p_ref               [1 -1 -2 0 0 0 0] 101325;  // Pa

    // Dimensionless parameters (no brackets):
    evaporationCoeff    0.05;
    alphaMin           0.001;

    // More dimensional parameters:
    maxSource          [1 -1 -3 0 0 0 0] 1e25;     // W/m³
}
```

## Quick Reference: Common Dimensions

| Physical Quantity | Dimensions | Example |
|------------------|------------|---------|
| Pressure (Pa) | `[1 -1 -2 0 0 0 0]` | `p_ref [1 -1 -2 0 0 0 0] 101325;` |
| Specific Energy (J/kg) | `[0 2 -2 0 0 0 0]` | `hf [0 2 -2 0 0 0 0] 9.1e6;` |
| Specific Gas Constant (J/(kg·K)) | `[0 2 -2 -1 0 0 0]` | `gasConstant [0 2 -2 -1 0 0 0] 174;` |
| Temperature (K) | `[0 0 0 1 0 0 0]` | `Tvap [0 0 0 1 0 0 0] 3560;` |
| Power Density (W/m³) | `[1 -1 -3 0 0 0 0]` | `maxSource [1 -1 -3 0 0 0 0] 1e25;` |
| Time (s) | `[0 0 1 0 0 0 0]` | `dt [0 0 1 0 0 0 0] 1e-15;` |
| Heat Capacity (J/(m³·K)) | `[1 -1 -2 -1 0 0 0]` | `Cl [1 -1 -2 -1 0 0 0] 2.5e6;` |
| Thermal Diffusivity (m²/s) | `[0 2 -1 0 0 0 0]` | `De [0 2 -1 0 0 0 0] 1e-4;` |

## Dimension Array Meaning

The seven numbers in brackets represent:
```
[kg  m  s  K  mol  A  cd]
 │   │  │  │   │   │  └─ Luminous intensity (candela) - rarely used
 │   │  │  │   │   └──── Electric current (ampere)
 │   │  │  │   └──────── Amount of substance (mole)
 │   │  │  └──────────── Temperature (kelvin)
 │   │  └─────────────── Time (second)
 │   └────────────────── Length (meter)
 └────────────────────── Mass (kilogram)
```

### Examples:
- Pressure: Pa = kg/(m·s²) = `[1 -1 -2 0 0 0 0]`
- Energy per mass: J/kg = m²/s² = `[0 2 -2 0 0 0 0]`
- Specific gas constant: J/(kg·K) = m²/(s²·K) = `[0 2 -2 -1 0 0 0]`

## How to Fix Your LIFT2 Case

1. Navigate to your case directory:
   ```bash
   cd ~/OpenFOAM/mavadi-v2406/run/LIFT2
   ```

2. Edit the controlDict file:
   ```bash
   gedit system/controlDict &
   # or
   nano system/controlDict
   ```

3. Find the `phaseChangeCoeffs` section (around lines 60-78)

4. Replace the entire section with the corrected version from:
   `/home/user/compInterFoam/FIXED_phaseChangeCoeffs_reference.txt`

5. Save and run again:
   ```bash
   compInterFoam
   ```

## Validation Checklist

Before running your simulation, check:

- [ ] All physical parameters have dimensional brackets
- [ ] No equals signs (=) in parameter definitions
- [ ] No colons (:) after parameter names
- [ ] Semicolons (;) at the end of each line
- [ ] Correct bracket format: `[# # # # # # #]`
- [ ] Dimensionless parameters (like coefficients) have NO brackets

## Working Examples

See the working test cases in this repository:
- `/home/user/compInterFoam/TEST1/system/controlDict` (lines 58-79)
- `/home/user/compInterFoam/TEST2/system/controlDict` (lines 58-80)
- `/home/user/compInterFoam/TestCase/system/controlDict`

## Need More Help?

If you encounter other errors:
1. Check the error message for the expected dimensions
2. Refer to the dimension reference table above
3. Compare with the working test cases
4. Verify the syntax matches exactly (no extra characters)
