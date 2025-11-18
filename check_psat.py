#!/usr/bin/env python3
"""
Titanium Vapor Pressure and Mass Flux Calculator
Verifies Clausius-Clapeyron thermodynamics for fs-LIFT simulation
"""
import math

# From controlDict phaseChangeCoeffs
L = 9.1e6      # J/kg (latent heat of vaporization for Ti)
R = 174        # J/kg/K (specific gas constant for Ti vapor)
T_vap = 2200   # K (reference temperature - CRITICAL PARAMETER)
p_ref = 101325 # Pa (atmospheric pressure - reference)
evap_coeff = 0.03  # Titanium evaporation coefficient

# Test temperatures from simulation
temperatures = [300, 1941, 2200, 3000, 5000, 7000, 9000, 9500]

print("=" * 80)
print("TITANIUM VAPOR PRESSURE & MASS FLUX DIAGNOSTIC")
print("=" * 80)
print(f"\nInput Parameters:")
print(f"  L (latent heat)     = {L:.2e} J/kg")
print(f"  R (gas constant)    = {R} J/kg/K")
print(f"  T_vap (reference)   = {T_vap} K  ← CRITICAL")
print(f"  p_ref (reference)   = {p_ref} Pa")
print(f"  evaporationCoeff    = {evap_coeff}")
print("\n" + "=" * 80)
print(f"{'T [K]':>8} {'Exponent':>12} {'p_sat [Pa]':>15} {'j_evap [kg/m²/s]':>18} {'Status':>20}")
print("=" * 80)

# Thresholds
massRateEps_old = 1e-8
massRateEps_new = 1e-12

sqrt_2pi = math.sqrt(2.0 * math.pi)

for T in temperatures:
    if T < 300:
        continue

    exponent = (L/R) * (1.0/T_vap - 1.0/T)
    p_sat = p_ref * math.exp(exponent)

    sqrt_R = math.sqrt(R)
    sqrt_T = math.sqrt(T)

    # Hertz-Knudsen equation
    j_evap = evap_coeff * p_sat / (sqrt_2pi * sqrt_R * sqrt_T)

    # Determine status
    if j_evap < massRateEps_new:
        status = "FILTERED (1e-12)"
    elif j_evap < massRateEps_old:
        status = "Was filtered (1e-8)"
    else:
        status = "ACTIVE ✓"

    print(f"{T:8.0f} {exponent:12.2f} {p_sat:15.2e} {j_evap:18.2e} {status:>20}")

print("=" * 80)
print("\nANALYSIS:")
print("-" * 80)

# Check at 9000K specifically
T_hot = 9000
exponent_hot = (L/R) * (1.0/T_vap - 1.0/T_hot)
p_sat_hot = p_ref * math.exp(exponent_hot)
j_evap_hot = evap_coeff * p_sat_hot / (sqrt_2pi * math.sqrt(R) * math.sqrt(T_hot))

print(f"\nAt T = 9000 K (your simulation max):")
print(f"  p_sat    = {p_sat_hot:.2e} Pa")
print(f"  j_evap   = {j_evap_hot:.2e} kg/m²/s")

if j_evap_hot < 1e-12:
    print(f"\n  *** PROBLEM IDENTIFIED ***")
    print(f"  Mass flux is BELOW threshold (1e-12) even at 9000K!")
    print(f"  This explains zero mass transfer.")
    print(f"\n  ROOT CAUSE: T_vap = {T_vap} K is too low for p_ref = 101325 Pa")
    print(f"\n  SOLUTION: Use T_vap = 3560 K (Ti boiling point at 1 atm)")
    print(f"            OR adjust p_ref to match T_vap = 2200 K")
    print(f"\n  Ti vapor pressure at 2200K is approximately: ~10-100 Pa (not 1 atm!)")

    # Calculate what p_ref should be
    # At boiling point (3560K), p_sat should equal p_ref
    print(f"\n  RECOMMENDED FIX:")
    print(f"  Option 1: Change T_vap to 3560 K (boiling point)")
    print(f"  Option 2: Change p_ref to ~10 Pa (matches T_vap=2200K)")

elif j_evap_hot < 1e-8:
    print(f"\n  Mass flux is detectable with new threshold (1e-12)")
    print(f"  Your fix (lowering massRateEps) should work!")

else:
    print(f"\n  Mass flux is HUGE - thermodynamics are correct!")
    print(f"  Problem must be elsewhere (activation windows, time gating, etc.)")

# Literature comparison
print(f"\n" + "-" * 80)
print(f"LITERATURE COMPARISON:")
print(f"-" * 80)
print(f"  Code comment says: 'Expected at 5000K: ~90 kg/m²/s'")
T_lit = 5000
exponent_lit = (L/R) * (1.0/T_vap - 1.0/T_lit)
p_sat_lit = p_ref * math.exp(exponent_lit)
j_lit = evap_coeff * p_sat_lit / (sqrt_2pi * math.sqrt(R) * math.sqrt(T_lit))
print(f"  Our calculation at 5000K: {j_lit:.2e} kg/m²/s")

if abs(j_lit - 90) / 90 < 0.5:
    print(f"  ✓ Within 50% of expected value - parameters reasonable")
elif j_lit < 0.01:
    print(f"  ✗ OFF BY ORDERS OF MAGNITUDE - parameters are wrong!")
    print(f"  ✗ Check T_vap and p_ref relationship")
else:
    print(f"  ~ Different but same order of magnitude")

# Physical validation
print(f"\n" + "=" * 80)
print(f"PHYSICAL VALIDATION:")
print(f"=" * 80)
print(f"\nTitanium Physical Properties (Literature):")
print(f"  Melting point:  1941 K")
print(f"  Boiling point:  3560 K (at 1 atm)")
print(f"  At 2200K: Superheated liquid (metastable)")
print(f"           Vapor pressure ~ 10-100 Pa (NOT 1 atm)")
print(f"\nClausius-Clapeyron Reference States:")
print(f"  Correct:   (T_vap=3560K, p_ref=101325 Pa) - boiling point")
print(f"  Incorrect: (T_vap=2200K, p_ref=101325 Pa) - mismatched!")
print(f"\n  Your current config uses mismatched reference state!")
print(f"=" * 80)

print(f"\nCONCLUSION:")
if j_evap_hot < 1e-10:
    print(f"  Thermodynamic parameters are WRONG.")
    print(f"  Must fix T_vap/p_ref relationship in controlDict.")
elif j_evap_hot < 1e-8:
    print(f"  Thermodynamics OK, threshold was too high.")
    print(f"  Your massRateEps fix (1e-12) should work.")
else:
    print(f"  Thermodynamics are CORRECT.")
    print(f"  Problem is elsewhere - check OpenFOAM code implementation.")
print(f"=" * 80)
