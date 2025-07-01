# Axial Compressor Design and Analysis Project

## Description

This MATLAB project designs and analyzes a multi-stage axial compressor for gas turbine applications. The project includes:
- Thermodynamic analysis of compressor stages
- Velocity triangle calculations
- Efficiency optimization
- Blade geometry design
- Performance parameter calculations

## Key Features

- **Efficiency Optimization**:
  - Iterative search for optimal flow coefficients (φ), work coefficients (ψ), and reaction degrees (R)
  - Calculates stage efficiency using thermodynamic relationships

- **Velocity Triangle Analysis**:
  - Computes inlet/outlet angles for rotors and stators
  - Determines absolute and relative velocities
  - Calculates Mach numbers throughout stages

- **Thermodynamic Calculations**:
  - Temperature and pressure ratios across stages
  - Total and static property calculations
  - Work input calculations

- **Geometric Design**:
  - Blade count determination
  - Chord length calculations
  - Compressor length estimation
  - Blade angle calculations (incidence, deviation, camber)

## Getting Started

### Prerequisites

- MATLAB (tested on R2019b and later)
- Symbolic Math Toolbox (for blade angle calculations)

### Installation

1. Place both files (`CompressorDesign_11.m` and `efficiency.m`) in the same directory
2. Run `CompressorDesign_11.m`

## Usage

The main script automatically performs all calculations and displays key results:

1. **Optimization Phase**:
   - Searches for optimal φ, ψ, R combination
   - Stores maximum efficiency configuration

2. **Performance Calculations**:
   - Velocity triangles
   - Temperature and pressure changes
   - Work input and efficiency

3. **Geometric Design**:
   - Number of blades
   - Compressor length
   - Blade angles

### Key Outputs

- Optimal efficiency and corresponding parameters
- Number of stages required
- Total power consumption
- Overall compressor efficiency
- Blade angles and geometry

## File Descriptions

- `CompressorDesign_11.m`: Main design and analysis script
- `efficiency.m`: Function calculating stage efficiency for given parameters

## Customization Options

### Input Parameters

Modify these in `CompressorDesign_11.m`:
```matlab
pressure_ratio = 3;    % Design pressure ratio
m_dot = 0.5;          % Mass flow rate (kg/s)
N = 30000;            % Rotational speed (RPM)
M_rel_max = 0.75;     % Maximum relative Mach number
DF_max = 0.5;         % Maximum diffusion factor
zeta = 0.75;          % Hub-to-tip ratio
```

### Optimization Ranges

Adjust parameter search spaces:
```matlab
phi_range = linspace(0.4,0.6,10);  % Flow coefficient
psi_range = linspace(0.2,0.4,10);  % Work coefficient
R_range = linspace(0.5,0.8,10);    % Reaction degree
```

## Technical Background

### Key Equations

1. **Stage Efficiency**:
   ```
   η = [(Pr^((γ-1)/γ)) - 1] / (τ - 1)
   ```
   Where Pr is pressure ratio and τ is temperature ratio

2. **Velocity Triangles**:
   - Calculated from flow coefficients and reaction degree
   - Used to determine work input and Mach numbers

3. **Blade Geometry**:
   - Uses Carter's rule for deviation angle prediction
   - Calculates camber and stagger angles

## Results Interpretation

- **Efficiency**: Higher values indicate better thermodynamic performance
- **Stage Count**: More stages may be needed for higher pressure ratios
- **Blade Angles**: Affect flow turning and diffusion
- **Mach Numbers**: Should remain below critical values to avoid losses

## License

This project is available for academic and research purposes. For commercial use, please contact the author.

## Acknowledgments

This implementation demonstrates fundamental gas turbine compressor design principles based on axial flow compressor theory. The methodology follows standard turbomachinery design practices.
