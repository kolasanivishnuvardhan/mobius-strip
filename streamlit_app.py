import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import integrate
import io
from app import MobiusStrip

# Set page configuration
st.set_page_config(
    page_title="Mobius Strip Explorer",
    page_icon="🔄",
    layout="wide"
)

# Title and description
st.title("🔄 Mobius Strip Explorer")
st.markdown("""
This interactive app allows you to explore the mathematical properties of a Mobius strip.
Adjust the parameters and see how they affect the geometry of this fascinating non-orientable surface.
""")

# Create sidebar for parameters
st.sidebar.header("Parameters")
radius = st.sidebar.slider("Radius (R)", 1.0, 5.0, 2.0, 0.1)
width = st.sidebar.slider("Width (w)", 0.2, 2.0, 1.0, 0.1)
resolution = st.sidebar.slider("Resolution", 20, 100, 50, 5)

# Visualization options
st.sidebar.header("Visualization Options")
show_wireframe = st.sidebar.checkbox("Show Wireframe", True)
show_edge = st.sidebar.checkbox("Highlight Edge", True)
alpha = st.sidebar.slider("Surface Transparency", 0.1, 1.0, 0.7, 0.1)

# Create tabs
tab1, tab2, tab3 = st.tabs(["3D Visualization", "Properties", "About"])

# Create Mobius strip instance
mobius = MobiusStrip(radius=radius, width=width, resolution=resolution)

# Compute properties
with st.spinner("Computing surface area..."):
    surface_area = mobius.compute_surface_area()
with st.spinner("Computing edge length..."):
    edge_length = mobius.compute_edge_length()

# Tab 1: 3D Visualization
with tab1:
    # Create matplotlib figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot the surface
    X, Y, Z = mobius.get_mesh_points()
    surface = ax.plot_surface(
        X, Y, Z,
        alpha=alpha,
        cmap='viridis',
        linewidth=0.1 if show_wireframe else 0,
        edgecolor='gray' if show_wireframe else 'none'
    )
    
    # Plot the edge if requested
    if show_edge:
        x_edge, y_edge, z_edge = mobius.get_edge_points()
        ax.plot(x_edge, y_edge, z_edge, 'r-', linewidth=3, label='Edge')
    
    # Set labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Mobius Strip (R={radius}, w={width})')
    
    # Equal aspect ratio
    max_range = max(
        np.ptp(X), np.ptp(Y), np.ptp(Z)
    ) / 2
    mid_x = (np.max(X) + np.min(X)) / 2
    mid_y = (np.max(Y) + np.min(Y)) / 2
    mid_z = (np.max(Z) + np.min(Z)) / 2
    
    ax.set_xlim(mid_x - max_range, mid_x + max_range)
    ax.set_ylim(mid_y - max_range, mid_y + max_range)
    ax.set_zlim(mid_z - max_range, mid_z + max_range)
    
    # Display the plot in Streamlit
    st.pyplot(fig)
    
    # Add rotation controls
    st.caption("Note: Interactive 3D rotation is limited in Streamlit. For full interactivity, consider downloading the Python script.")
    
    # Add download button for the figure
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    st.download_button(
        label="Download 3D Visualization",
        data=buf,
        file_name="mobius_strip.png",
        mime="image/png"
    )

# Tab 2: Properties
with tab2:
    # Create two columns
    col1, col2 = st.columns(2)
    
    # Column 1: Display computed properties
    with col1:
        st.subheader("Computed Properties")
        st.metric("Surface Area", f"{surface_area:.4f} square units")
        st.metric("Edge Length", f"{edge_length:.4f} units")
        
        st.subheader("Parameters")
        st.write(f"Radius (R): {radius}")
        st.write(f"Width (w): {width}")
        st.write(f"Resolution: {resolution}")
        st.write(f"Parameter ranges:")
        st.write(f"- u ∈ [0, 2π] ≈ [0, {2*np.pi:.3f}]")
        st.write(f"- v ∈ [-w/2, w/2] = [{-width/2:.3f}, {width/2:.3f}]")
    
    # Column 2: Show formulas and mathematical background
    with col2:
        st.subheader("Mathematical Model")
        st.markdown(r"""
        The Mobius strip is defined by the parametric equations:
        
        $$x(u,v) = (R + v \cdot \cos(u/2)) \cdot \cos(u)$$
        
        $$y(u,v) = (R + v \cdot \cos(u/2)) \cdot \sin(u)$$
        
        $$z(u,v) = v \cdot \sin(u/2)$$
        
        Where:
        - $u \in [0, 2\pi]$ (parameter along the strip)
        - $v \in [-w/2, w/2]$ (parameter across the width)
        
        The surface area is computed using the double integral:
        
        $$A = \iint \left\| \frac{\partial \vec{r}}{\partial u} \times \frac{\partial \vec{r}}{\partial v} \right\| du dv$$
        
        The edge length is computed by integrating along the boundary:
        
        $$L = \int_0^{2\pi} \left\| \frac{d\vec{r}}{du} \right\| du$$
        
        where $\vec{r}$ is evaluated at $v = w/2$ (the edge of the strip).
        """)

# Tab 3: About
with tab3:
    st.subheader("About the Mobius Strip")
    st.markdown("""
    The Mobius strip (or Möbius band) is a surface with only one side and one edge. It can be created by taking a strip of paper, giving it a half-twist, and joining the ends.
    
    ### Key Properties:
    - **Non-orientable surface**: It has only one side
    - **Single edge**: Despite appearing to have two edges, it has only one continuous edge
    - **Fascinating topology**: Cutting a Mobius strip along the center line results in a single longer strip with two full twists
    
    ### Applications:
    - Mathematical topology and research
    - Engineering (conveyor belts that wear evenly on both sides)
    - Art and design
    - Puzzles and recreational mathematics
    """)
    
    st.subheader("About this App")
    st.markdown("""
    This app was created using:
    - Python for mathematical modeling
    - NumPy and SciPy for numerical calculations
    - Matplotlib for 3D visualization
    - Streamlit for the web interface
    
    The source code for this app is available on GitHub.
    """)

# Add a footer
st.markdown("---")
st.markdown("Created with ❤️ using Streamlit and Python") 