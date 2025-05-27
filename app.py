import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate
import warnings
warnings.filterwarnings('ignore')

class MobiusStrip:
    """
    A class to model and analyze a Mobius strip using parametric equations.

    The Mobius strip is defined by the parametric equations:
    x(u,v) = (R + v*cos(u/2)) * cos(u)
    y(u,v) = (R + v*cos(u/2)) * sin(u)
    z(u,v) = v * sin(u/2)

    Where:
    - u ∈ [0, 2π] (parameter along the strip)
    - v ∈ [-w/2, w/2] (parameter across the width)
    """

    def __init__(self, radius=2.0, width=1.0, resolution=50):
        """
        Initialize the Mobius strip with given parameters.

        Parameters:
        -----------
        radius : float
            Distance from center to the strip (R)
        width : float
            Width of the strip (w)
        resolution : int
            Number of points in each parameter direction for mesh generation
        """
        self.radius = radius
        self.width = width
        self.resolution = resolution

        # Generate parameter grids
        self.u = np.linspace(0, 2*np.pi, resolution)
        self.v = np.linspace(-width/2, width/2, resolution)
        self.U, self.V = np.meshgrid(self.u, self.v)

        # Compute the 3D mesh
        self._compute_mesh()

    def _compute_mesh(self):
        """Compute the 3D coordinates of the Mobius strip surface."""
        # Parametric equations for Mobius strip
        self.X = (self.radius + self.V * np.cos(self.U/2)) * np.cos(self.U)
        self.Y = (self.radius + self.V * np.cos(self.U/2)) * np.sin(self.U)
        self.Z = self.V * np.sin(self.U/2)

    def get_mesh_points(self):
        """
        Return the mesh points as arrays.

        Returns:
        --------
        tuple: (X, Y, Z) arrays representing the 3D mesh
        """
        return self.X, self.Y, self.Z

    def parametric_point(self, u, v):
        """
        Calculate a single point on the Mobius strip for given parameters.

        Parameters:
        -----------
        u : float or array
            Parameter along the strip [0, 2π]
        v : float or array
            Parameter across the width [-w/2, w/2]

        Returns:
        --------
        tuple: (x, y, z) coordinates
        """
        x = (self.radius + v * np.cos(u/2)) * np.cos(u)
        y = (self.radius + v * np.cos(u/2)) * np.sin(u)
        z = v * np.sin(u/2)
        return x, y, z

    def _partial_derivatives(self, u, v):
        """
        Calculate partial derivatives for surface area computation.

        Returns:
        --------
        tuple: (du_vec, dv_vec) - partial derivative vectors
        """
        # Partial derivatives with respect to u
        dx_du = -(self.radius + v * np.cos(u/2)) * np.sin(u) - (v/2) * np.sin(u/2) * np.cos(u)
        dy_du = (self.radius + v * np.cos(u/2)) * np.cos(u) - (v/2) * np.sin(u/2) * np.sin(u)
        dz_du = (v/2) * np.cos(u/2)

        # Partial derivatives with respect to v
        dx_dv = np.cos(u/2) * np.cos(u)
        dy_dv = np.cos(u/2) * np.sin(u)
        dz_dv = np.sin(u/2)

        du_vec = np.array([dx_du, dy_du, dz_du])
        dv_vec = np.array([dx_dv, dy_dv, dz_dv])

        return du_vec, dv_vec

    def _surface_element_magnitude(self, u, v):
        """
        Calculate the magnitude of the cross product of partial derivatives.
        This gives the differential surface area element.
        """
        du_vec, dv_vec = self._partial_derivatives(u, v)

        # Cross product
        cross_product = np.cross(du_vec.T, dv_vec.T, axis=-1)

        # Magnitude of cross product
        magnitude = np.sqrt(np.sum(cross_product**2, axis=-1))

        return magnitude

    def compute_surface_area(self):
        """
        Compute the surface area of the Mobius strip using numerical integration.

        The surface area is given by the double integral:
        ∬ ||∂r/∂u × ∂r/∂v|| du dv

        Returns:
        --------
        float: Approximate surface area
        """
        try:
            # Define the integrand function
            def integrand(v, u):
                return self._surface_element_magnitude(u, v)

            # Perform double integration
            result, error = integrate.dblquad(
                integrand,
                0, 2*np.pi,  # u limits
                lambda u: -self.width/2, lambda u: self.width/2,  # v limits
                epsabs=1e-6, epsrel=1e-6
            )

            self.surface_area = result
            self.surface_area_error = error

        except Exception as e:
            print(f"Integration failed: {e}")
            # Fallback to grid-based approximation
            self.surface_area = self._approximate_surface_area_grid()
            self.surface_area_error = None

        return self.surface_area

    def _approximate_surface_area_grid(self):
        """
        Approximate surface area using a discrete grid method.
        This is a fallback method when numerical integration fails.
        """
        du = 2*np.pi / (self.resolution - 1)
        dv = self.width / (self.resolution - 1)

        total_area = 0
        for i in range(self.resolution - 1):
            for j in range(self.resolution - 1):
                u_val = self.u[i]
                v_val = self.v[j]

                # Calculate surface element at this point
                magnitude = self._surface_element_magnitude(u_val, v_val)
                total_area += magnitude * du * dv

        return total_area

    def compute_edge_length(self):
        """
        Compute the length of the edge of the Mobius strip.

        The edge consists of two parts:
        1. Outer edge: v = w/2
        2. Inner edge: v = -w/2 (but they form a single continuous edge)

        Returns:
        --------
        float: Total edge length
        """
        def edge_speed(u, v_fixed):
            """Calculate the speed along the edge parameterized by u."""
            du_vec, _ = self._partial_derivatives(u, v_fixed)
            return np.sqrt(np.sum(du_vec**2))

        # Integrate along the edge (the Mobius strip has only one edge!)
        edge_length, _ = integrate.quad(
            lambda u: edge_speed(u, self.width/2),
            0, 2*np.pi,
            epsabs=1e-6, epsrel=1e-6
        )

        self.edge_length = edge_length
        return edge_length

    def get_edge_points(self, num_points=100):
        """
        Get points along the edge of the Mobius strip for visualization.

        Parameters:
        -----------
        num_points : int
            Number of points to generate along the edge

        Returns:
        --------
        tuple: (x, y, z) coordinates of edge points
        """
        u_edge = np.linspace(0, 2*np.pi, num_points)
        v_edge = self.width/2  # Outer edge

        x_edge, y_edge, z_edge = self.parametric_point(u_edge, v_edge)

        return x_edge, y_edge, z_edge

    def visualize(self, show_wireframe=True, show_edge=True, alpha=0.7):
        """
        Create a 3D visualization of the Mobius strip.

        Parameters:
        -----------
        show_wireframe : bool
            Whether to show wireframe
        show_edge : bool
            Whether to highlight the edge
        alpha : float
            Transparency of the surface
        """
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Plot the surface
        surface = ax.plot_surface(
            self.X, self.Y, self.Z,
            alpha=alpha,
            cmap='viridis',
            linewidth=0.1 if show_wireframe else 0,
            edgecolor='gray' if show_wireframe else 'none'
        )

        # Plot the edge if requested
        if show_edge:
            x_edge, y_edge, z_edge = self.get_edge_points()
            ax.plot(x_edge, y_edge, z_edge, 'r-', linewidth=3, label='Edge')

        # Set labels and title
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title(f'Mobius Strip (R={self.radius}, w={self.width})')

        # Add colorbar
        plt.colorbar(surface, ax=ax, shrink=0.8)

        if show_edge:
            ax.legend()

        # Equal aspect ratio
        max_range = max(
            np.ptp(self.X), np.ptp(self.Y), np.ptp(self.Z)
        ) / 2
        mid_x = (np.max(self.X) + np.min(self.X)) / 2
        mid_y = (np.max(self.Y) + np.min(self.Y)) / 2
        mid_z = (np.max(self.Z) + np.min(self.Z)) / 2

        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)

        plt.tight_layout()
        plt.show()

    def print_properties(self):
        """Print computed geometric properties."""
        print("="*50)
        print("MOBIUS STRIP GEOMETRIC PROPERTIES")
        print("="*50)
        print(f"Parameters:")
        print(f"  Radius (R): {self.radius}")
        print(f"  Width (w): {self.width}")
        print(f"  Resolution: {self.resolution}")
        print(f"\nComputed Properties:")

        if hasattr(self, 'surface_area'):
            print(f"  Surface Area: {self.surface_area:.6f}")
            if self.surface_area_error:
                print(f"  Integration Error: ±{self.surface_area_error:.2e}")

        if hasattr(self, 'edge_length'):
            print(f"  Edge Length: {self.edge_length:.6f}")

        print(f"\nMesh Information:")
        print(f"  Total mesh points: {self.resolution}×{self.resolution} = {self.resolution**2}")
        print(f"  Parameter ranges:")
        print(f"    u ∈ [0, 2π] ≈ [0, {2*np.pi:.3f}]")
        print(f"    v ∈ [-w/2, w/2] = [{-self.width/2:.3f}, {self.width/2:.3f}]")


def demo_mobius_analysis():
    """
    Demonstration function showing various Mobius strip analyses.
    """
    print("MOBIUS STRIP MATHEMATICAL MODELING DEMO")
    print("=" * 60)

    # Create different Mobius strips for comparison
    configs = [
        {"radius": 2.0, "width": 1.0, "resolution": 50, "name": "Standard"},
        {"radius": 3.0, "width": 0.5, "resolution": 60, "name": "Narrow"},
        {"radius": 1.5, "width": 1.5, "resolution": 40, "name": "Wide"}
    ]

    results = []

    for config in configs:
        print(f"\nAnalyzing {config['name']} Mobius Strip...")

        # Create Mobius strip
        mobius = MobiusStrip(
            radius=config["radius"],
            width=config["width"],
            resolution=config["resolution"]
        )

        # Compute properties
        surface_area = mobius.compute_surface_area()
        edge_length = mobius.compute_edge_length()

        results.append({
            "name": config["name"],
            "radius": config["radius"],
            "width": config["width"],
            "surface_area": surface_area,
            "edge_length": edge_length,
            "mobius": mobius
        })

    # Print comparison table
    print("\n" + "="*80)
    print("COMPARISON OF DIFFERENT MOBIUS STRIPS")
    print("="*80)
    print(f"{'Configuration':<12} {'Radius':<8} {'Width':<8} {'Surf. Area':<12} {'Edge Length':<12}")
    print("-" * 80)

    for result in results:
        print(f"{result['name']:<12} {result['radius']:<8.1f} {result['width']:<8.1f} "
              f"{result['surface_area']:<12.4f} {result['edge_length']:<12.4f}")

    # Visualize the first configuration
    print(f"\nVisualizing {results[0]['name']} Mobius Strip...")
    results[0]['mobius'].print_properties()
    results[0]['mobius'].visualize()

    return results


if __name__ == "__main__":
    # Run the demonstration
    results = demo_mobius_analysis()

    # Additional analysis: Show how properties scale with parameters
    print("\n" + "="*60)
    print("PARAMETER SCALING ANALYSIS")
    print("="*60)

    # Test different radii with fixed width
    print("\nSurface area scaling with radius (width = 1.0):")
    radii = [1.0, 1.5, 2.0, 2.5, 3.0]
    for r in radii:
        m = MobiusStrip(radius=r, width=1.0, resolution=30)
        area = m.compute_surface_area()
        print(f"  R = {r:.1f}: Area = {area:.4f}")

    # Test different widths with fixed radius
    print("\nSurface area scaling with width (radius = 2.0):")
    widths = [0.5, 0.75, 1.0, 1.25, 1.5]
    for w in widths:
        m = MobiusStrip(radius=2.0, width=w, resolution=30)
        area = m.compute_surface_area()
        print(f"  w = {w:.2f}: Area = {area:.4f}")