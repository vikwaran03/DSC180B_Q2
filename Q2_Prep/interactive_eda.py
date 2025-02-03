import numpy as np
import pandas as pd
import networkx as nx
from pyvis.network import Network
import streamlit as st

# Set Streamlit page to wide layout (must be at the very start)
st.set_page_config(layout="wide")

# Add at the beginning after st.set_page_config
if 'selected_node' not in st.session_state:
    st.session_state.selected_node = None

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
        max-width: 100% !important; /* Ensures full width */
        padding: 0 !important;
        margin: 0 !important;
        background-color: white !important;
        color: black !important; /* Ensures black font */
    }
    iframe {
        width: 100%;
        height: 90vh;
        margin: 0;
        padding: 0;
        border: none; /* Removes iframe border */
    }
    header, footer {
        display: none !important;
    }
    /* Ensures all text elements are black */
    h1, h2, h3, h4, h5, h6, p, span, div, label, .st-emotion-cache, .stMarkdown, .stTextInput, .stButton, .stSelectbox, .stSlider, .stCheckbox, .stRadio {
        color: black !important;
    }
    /* Ensures Streamlit widgets have black text */
    .st-bb, .st-at, .st-ae, .st-af, .st-ag, .st-ah, .st-ai, .st-aj, .st-ak, .st-al, .st-am, .st-an, .st-ao, .st-ap, .st-aq, .st-ar, .st-as {
        color: black !important;
    }
    /* Ensures no extra padding or margins in Streamlit containers */
    .st-emotion-cache-1v0mbdj, .st-emotion-cache-1y4p8pa, .st-emotion-cache-1n76uvr {
        padding: 0 !important;
        margin: 0 !important;
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
    
toggle_script = """
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function () {
    initializeNetwork();
});

function initializeNetwork() {
    storeOriginalColors();
    setupNetworkListeners();
    restoreSelection(true); // Apply selection *before* any network updates
}

function preventResetOnUpdate() {
    const selectedNode = sessionStorage.getItem("selectedNode");
    if (selectedNode !== null) {
        // Apply selection immediately before any visual updates
        requestAnimationFrame(() => {
            toggleEdgesAndNodes(parseInt(selectedNode), true);
        });
    }
}

function toggleEdgesAndNodes(selectedNode, preventSave = false) {
    const network = window.network;
    if (!network) return;
    
    // Store positions before update
    network.storePositions();
    
    const nodes = network.body.data.nodes;
    const edges = network.body.data.edges;
    
    if (selectedNode === null) {
        nodes.update(nodes.get().map(node => ({
            id: node.id,
            color: node.originalColor
        })));
        edges.update(edges.get().map(edge => ({ ...edge, hidden: false })));
    } else {
        const connectedNodes = network.getConnectedNodes(selectedNode);
        const connectedEdges = network.getConnectedEdges(selectedNode);
        
        // Update nodes immediately with new state
        nodes.update(nodes.get().map(node => ({
            id: node.id,
            color: (node.id === selectedNode || connectedNodes.includes(node.id)) 
                ? node.originalColor 
                : { background: '#CCCCCC', border: '#CCCCCC' }
        })));
        
        edges.update(edges.get().map(edge => ({
            ...edge, 
            hidden: !connectedEdges.includes(edge.id)
        })));
        
        if (!preventSave) {
            sessionStorage.setItem("selectedNode", selectedNode);
        }
    }
}


// Restore selection instantly before updates occur
function restoreSelection(preventSave = false) {
    var selectedNode = sessionStorage.getItem("selectedNode");
    if (selectedNode !== null && window.network) {
        toggleEdgesAndNodes(parseInt(selectedNode), preventSave);
    }
}

// Ensure original colors persist after network loads
function storeOriginalColors() {
    if (window.network && window.network.body && window.network.body.data) {
        var nodes = window.network.body.data.nodes;
        nodes.update(nodes.get().map(node => ({
            id: node.id,
            originalColor: node.color || { background: '#97C2FC', border: '#2B7CE9' }
        })));
    }
}

// Setup event listeners
function setupNetworkListeners() {
    var network = window.network;
    if (!network) return;

    network.on("click", function(params) {
        if (params.nodes.length > 0) {
            toggleEdgesAndNodes(params.nodes[0]);
        } else {
            toggleEdgesAndNodes(null);
            sessionStorage.removeItem("selectedNode");
        }
    });

    // Hook into the slider update to prevent flickering
    var slider = document.getElementById("your-slider-id");
    if (slider) {
        slider.addEventListener("input", function () {
            preventResetOnUpdate();
        });
    }
}

// Prevent selection reset by applying the selection before the UI updates
function preventResetOnUpdate() {
    var selectedNode = sessionStorage.getItem("selectedNode");
    if (selectedNode !== null) {
        toggleEdgesAndNodes(parseInt(selectedNode), true);
    }
}
</script>
"""

zoom_control_script = """
<script type="text/javascript">
document.addEventListener("DOMContentLoaded", function() {
    let isAnimating = false;
    let initialScale, initialPosition, minScale;

    const setupZoom = () => {
        if (!window.network || typeof window.network.getScale !== "function") return;

        initialScale = window.network.getScale();
        initialPosition = window.network.getViewPosition();
        minScale = initialScale * 1;

        // Clear existing listeners to prevent duplicates
        window.network.off("zoom");
        window.network.off("dragEnd");
    };

    function resetView() {
        if (isAnimating) return;
        isAnimating = true;
        
        window.network.moveTo({
            scale: minScale,
            position: initialPosition,
            animation: {
                duration: 1000,
                easingFunction: "easeInOutCubic"
            }
        });
        
        setTimeout(() => {
            isAnimating = false;
        }, 1000);
    }

    const checkScaleAndReset = () => {
        const currentScale = window.network.getScale();
        if (currentScale <= minScale && !isAnimating) {
            resetView();
        }
    };

    const initInterval = setInterval(() => {
        if (window.network && typeof window.network.getScale === "function") {
            setupZoom();

            window.network.on("zoom", checkScaleAndReset);
            window.network.on("dragEnd", checkScaleAndReset);

            window.addEventListener('resize', setupZoom);

            clearInterval(initInterval);
        }
    }, 300);
});
</script>
"""


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
    source_code + toggle_script + zoom_control_script,
    height=1000,  # Increase the height here as needed
    scrolling=True
)