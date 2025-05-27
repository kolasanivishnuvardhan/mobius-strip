# Mobius Strip Mathematical Modeling

This project implements a mathematical model of a Mobius strip using parametric equations and computes its key geometric properties.

## Overview

The Mobius strip is a fascinating non-orientable surface with only one side and one edge. This implementation models the strip using parametric equations and provides tools to analyze its geometric properties.

## Features

- **Parametric Modeling**: Uses the standard parametric equations for the Mobius strip
- **Geometric Properties**: Computes surface area and edge length using numerical integration
- **Visualization**: 3D visualization with matplotlib including surface coloring and edge highlighting
- **Parameter Analysis**: Tools to analyze how properties scale with different parameters

## Mathematical Model

The Mobius strip is defined by the parametric equations:
- x(u,v) = (R + v·cos(u/2)) · cos(u)
- y(u,v) = (R + v·cos(u/2)) · sin(u)
- z(u,v) = v · sin(u/2)

Where:
- u ∈ [0, 2π] (parameter along the strip)
- v ∈ [-w/2, w/2] (parameter across the width)
- R = radius (distance from center to the strip)
- w = width of the strip

## Implementation Details

### Code Structure
- `MobiusStrip` class handles all calculations and visualization
- Numerical integration for surface area and edge length using scipy
- Fallback methods for approximation when integration fails
- Visualization tools with matplotlib

### Surface Area Approximation
The surface area is computed using the double integral:
∬ ||∂r/∂u × ∂r/∂v|| du dv

This is implemented using scipy's `dblquad` function, with a grid-based approximation as fallback.

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- SciPy

## Usage

```python
# Create a Mobius strip
mobius = MobiusStrip(radius=2.0, width=1.0, resolution=50)

# Compute properties
surface_area = mobius.compute_surface_area()
edge_length = mobius.compute_edge_length()

# Display properties
mobius.print_properties()

# Visualize
mobius.visualize()
```

## Challenges and Solutions

- **Numerical Integration**: The surface area calculation requires careful handling of the partial derivatives and cross product
- **Visualization**: Creating an effective 3D visualization with proper aspect ratio and coloring
- **Edge Representation**: The Mobius strip has only one continuous edge, which required special handling

## Sample Output

The code provides detailed analysis of different Mobius strip configurations, showing how geometric properties scale with parameters like radius and width. 