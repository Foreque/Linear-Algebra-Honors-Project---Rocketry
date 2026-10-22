# data/

Reference inputs for the simulator and generated datasets.

## Required inputs

| File | Used in | Source |
|---|---|---|
| `AeroTech_L1940X.eng` | Stage 4 (thrust curve, burn time, propellant mass) | ThrustCurve.org, AeroTech L1940X, RASP format, manufacturer data |
| `CD_Test.CSV` | Stage 5 (Cd vs. Mach lookup) | RasAero II export for Zenith 2 |

Motor reference values from ThrustCurve.org: total impulse 4329.2 N·s, average thrust 1939.3 N, burn time 2.2 s, total mass 3856 g, propellant mass 1825 g.

## Notes

- Record the RasAero II version and the export date when `CD_Test.CSV` is regenerated.
- The Cd curve that generates synthetic data is the ground truth for the twin experiment. Do not edit it by hand.
