
import numpy as np
import pandas as pd
import networkx as nx
from pyvis.network import Network
import streamlit as st

# Set Streamlit page to wide layout (must be at the very start)
st.set_page_config(layout="wide")

# Load data
ec_hic = np.load('data/GBM39ec_5k_collapsed_matrix.npy')
hsr_hic = np.load('data/GBM39HSR_5k_collapsed_matrix.npy')

hsr_df = pd.read_csv('data/HSR_features.csv')
hic_hsr_chr7 = hsr_df[(hsr_df['start'] >= 54765000) & (hsr_df['end'] <= 56050000)]
hic_hsr_chr7 = hic_hsr_chr7[hic_hsr_chr7['chromosome'] == 'NC_000007.14']

# Streamlit UI
st.title("Interactive Graph Visualization")

# Toggle option for node coloring
view_option = st.radio(
    "Select Node Attribute to Display",
    ["Total Genes", "Read Counts"],
    index=0
)

st.markdown(
    """
    <style>
    html, body, [class*="stApp"] {
        width: 100%;
        height: 100%;
        padding: 0;
        margin: 0;
        overflow: hidden;
        background-color: white !important; /* White background */
        color: black !important; /* Black text */
    }
    .main .block-container {
        width: 100%;
        padding: 0;
        margin: 0;
        background-color: white !important;
        color: black !important; /* Ensures black font */
    }
    iframe {
        width: 100%;
        height: 90vh;
        margin: 0;
        padding: 0;
        border: true;
    }
    header, footer {
        display: none !important;
    }
    .st-emotion-cache {
        color: black !important; /* Ensures black font for Streamlit elements */
    }
    h1, h2, h3, h4, h5, h6, p, span, div, label {
        color: black !important; /* Forces black font for all text */
    }
    </style>
    """,
    unsafe_allow_html=True
)


# Add a slider for threshold filtering
to_filter = st.slider(
    "Percentile Threshold (to_filter)",
    min_value=70,
    max_value=95,
    value=90,  # Default value
    step=5,
    help="Adjust the percentile threshold for filtering the HSR Hi-C matrix."
)

# Filter Hi-C matrix based on the slider value
thres = np.percentile(hsr_hic, to_filter)
vis_hsr = hsr_hic.copy()
vis_hsr[vis_hsr < thres] = 0
np.fill_diagonal(vis_hsr, 0)

# Create graph
G = nx.from_numpy_array(vis_hsr)

# Normalize edge weights
edge_weights = nx.get_edge_attributes(G, "weight")
max_weight = max(edge_weights.values())
min_weight = min(edge_weights.values())
normalized_weights = {
    edge: (weight - min_weight) * 0.3 / (max_weight - min_weight)
    for edge, weight in edge_weights.items()
}

# Node attributes
tot_genes = hic_hsr_chr7.total_genes.to_numpy()
num_reads = np.log(hic_hsr_chr7.read_count.to_numpy() + 0.0000000001)

# Create a Pyvis network with full-screen dimensions
net = Network(
    height="100vh",
    width="100%",
    bgcolor="#FFFFFF",
    font_color="black"
)

# Layout
pos = nx.circular_layout(G)

# Add nodes dynamically based on view selection
for node in G.nodes:
    x, y = pos[node]
    if view_option == "Total Genes":
        color = '#4B0082'
        if tot_genes[node] == 1:
            color = '#00CED1'
        elif tot_genes[node] >= 2:
            color = '#FFD700'
        title_text = f"Total Genes: {tot_genes[node]}"
    else:
        color = '#4B0082'
        if num_reads[node] > 1000:
            color = '#00CED1'
        elif num_reads[node] > 2000:
            color = '#FFD700'
        title_text = f"Read Counts: {num_reads[node]:.2f}"

    net.add_node(
        node,
        label=str(node),
        title=title_text,
        size=15,
        color={'background': color, 'border': color},
        x=x * 1000,
        y=y * 1000
    )

# Add edges
for edge in G.edges:
    net.add_edge(
        edge[0],
        edge[1],
        value=normalized_weights[edge],
        title=f"Weight: {edge_weights[edge]:.2f}",
        color="#D3D3D3",
        smooth=False
    )

# JavaScript for toggling edge visibility when clicking a node
toggle_script = """
<script type="text/javascript">
function toggleEdges(selectedNode) {
    var network = window.network;
    var edges = network.body.data.edges;
    
    if (selectedNode === null) {
        // Show all edges if no node is selected
        edges.update(edges.get().map(edge => ({ ...edge, hidden: false })));
    } else {
        var connectedEdges = network.getConnectedEdges(selectedNode);
        
        // Hide all edges first
        edges.update(edges.get().map(edge => ({ ...edge, hidden: true })));
        
        // Show only edges connected to the selected node
        edges.update(connectedEdges.map(edgeId => ({ id: edgeId, hidden: false })));
    }
}

// Listen for node selection
function setupNetworkListeners() {
    var network = window.network;
    network.on("click", function (params) {
        if (params.nodes.length > 0) {
            toggleEdges(params.nodes[0]);
        } else {
            toggleEdges(null);  // Restore all edges if clicking outside a node
        }
    });
}

// Initialize the listener when the page loads
document.addEventListener("DOMContentLoaded", setupNetworkListeners);
</script>
"""

# (Optional) Additional JavaScript for node highlight if needed (highlight_script omitted for brevity)

# Set options for network
net.set_options("""
{
    "nodes": {
        "borderWidth": 0,
        "shape": "dot",
        "font": {
            "size": 10,
            "color": "#000000",
            "face": "arial",
            "align": "center"
        }
    },
    "edges": {
        "width": 0.3,
        "selectionWidth": 0.3,
        "smooth": false,
        "color": {
            "opacity": 0.3
        }
    },
    "interaction": {
        "dragNodes": false,
        "dragView": true,
        "zoomView": true,
        "hover": true
    }, 
    "physics": {
        "enabled": false
    }
}
""")

# Save and display the graph
net.save_graph("graph.html")
HtmlFile = open("graph.html", "r", encoding="utf-8")
source_code = HtmlFile.read()

# Embed the graph with JavaScript and increase height parameter
st.components.v1.html(
    source_code + toggle_script,
    height=1000,  # Increase the height here as needed
    scrolling=True
)