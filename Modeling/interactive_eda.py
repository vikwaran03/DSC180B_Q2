import numpy as np
import pandas as pd
import networkx as nx
from pyvis.network import Network
import streamlit as st
import os

st.set_page_config(layout="wide")

if 'selected_node' not in st.session_state:
    st.session_state.selected_node = None

@st.cache_data
def get_data_path():
    return os.path.join(os.path.dirname(__file__), 'data')

@st.cache_data
def load_hic_data():
    data_path = get_data_path()
    ec_hic = np.load(os.path.join(data_path, 'GBM39ec_5k_collapsed_matrix.npy'))
    hsr_hic = np.load(os.path.join(data_path, 'GBM39HSR_5k_collapsed_matrix.npy'))
    return ec_hic, hsr_hic

@st.cache_data
def load_hsr_features():
    data_path = get_data_path()
    hsr_df = pd.read_csv(os.path.join(data_path, 'HSR_features.csv'))
    hic_hsr_chr7 = hsr_df[(hsr_df['start'] >= 54765000) & (hsr_df['end'] <= 56050000)]
    hic_hsr_chr7 = hic_hsr_chr7[hic_hsr_chr7['chromosome'] == 'NC_000007.14']
    return hic_hsr_chr7

@st.cache_data
def load_ec_features():
    data_path = get_data_path()
    ec_df = pd.read_csv(os.path.join(data_path, 'ecDNA_features.csv'))
    hic_ec_chr7 = ec_df[(ec_df['start'] >= 54765000) & (ec_df['end'] <= 56050000)]
    hic_ec_chr7 = hic_ec_chr7[hic_ec_chr7['chromosome'] == 'NC_000007.14']
    return hic_ec_chr7

data_path = get_data_path()
ec_hic, hsr_hic = load_hic_data()
hic_ec_chr7 = load_ec_features()
hic_hsr_chr7 = load_hsr_features()

st.title("Interactive Graph Visualization")

dataset_option = st.radio(
    "Select Dataset",
    ["HSR", "EC"],
    index=0,
    help="Choose between HSR or EC dataset for visualization."
)

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


if dataset_option == "HSR":
    selected_hic = hsr_hic
    selected_features = hic_hsr_chr7
else:
    selected_hic = ec_hic
    selected_features = hic_ec_chr7

to_filter = st.slider(
    "Percentile Threshold (to_filter)",
    min_value=70,
    max_value=95,
    value=90,  # Default value
    step=5,
)

thres = np.percentile(selected_hic, to_filter)
vis_matrix = selected_hic.copy()
vis_matrix[vis_matrix < thres] = 0
np.fill_diagonal(vis_matrix, 0)

@st.cache_data
def generate_graph(vis_matrix):
    G = nx.from_numpy_array(vis_matrix)
    return G

G = generate_graph(vis_matrix)

edge_weights = nx.get_edge_attributes(G, "weight")
if edge_weights:  
    max_weight = max(edge_weights.values())
    min_weight = min(edge_weights.values())
else:
    max_weight, min_weight = 1, 0 

normalized_weights = {
    edge: (weight - min_weight) * 0.3 / (max_weight - min_weight)
    for edge, weight in edge_weights.items()
}

tot_genes = selected_features.total_genes.to_numpy()
num_reads = np.log(selected_features.read_count.to_numpy() + 1e-10)

net = Network(
    height="100vh",
    width="100%",
    bgcolor="#FFFFFF",
    font_color="black"
)

pos = nx.circular_layout(G)

thresh1, thresh2 = np.percentile(num_reads, 50), np.percentile(num_reads, 75)

for node in G.nodes:
    tit_text = f"""
            Node: {node}
            {view_option}: {tot_genes[node] if view_option == 'Total Genes' else num_reads[node]:.2f}
            Connected Edges: {len(G.edges(node))}
            """
    
    x, y = pos[node]
    
    if view_option == "Total Genes":
        color = '#1f77b4'  # Base blue
        if tot_genes[node] == 1:
            color = '#2ca02c'  # Green
        elif tot_genes[node] >= 2:
            color = '#ff7f0e'  # Orange
        title_text = f"Total Genes: {tit_text}"
    
    else:  # View option is Read Counts
        color = '#1f77b4'
        if num_reads[node] > thresh2:
            color = '#ff7f0e'  # Orange
        elif num_reads[node] > thresh1:
            color = '#2ca02c'  # Green
        
        title_text = f"Read Counts: {tit_text}"

    net.add_node(
        node,
        label=str(node),
        title=title_text,
        size=15,
        color={'background': color, 'border': color},
        x=x * 1000,
        y=y * 1000,
    )

for edge in G.edges:
    net.add_edge(
        edge[0],
        edge[1],
        value=normalized_weights.get(edge, 0),
        title=f"Weight: {edge_weights.get(edge, 0):.2f}",
        color="#D3D3D3",
        smooth=False,
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
    let initialScale, initialPosition, minScale, maxScale;
    let isAnimating = false;
    // Define a threshold (in pixels) for how far the view can be dragged when near default zoom.
    const positionThreshold = 200; 

    const setupZoom = () => {
        if (!window.network || typeof window.network.getScale !== "function") return;

        // Record the current scale and position as baseline.
        initialScale = window.network.getScale();
        initialPosition = window.network.getViewPosition();
        minScale = initialScale * 1;      // Minimum zoom allowed (baseline)
        maxScale = initialScale * 1;      // Maximum zoom allowed (2Ã— baseline)

        // Enforce zoom bounds so the user cannot zoom in past maxScale.
        window.network.setOptions({
            interaction: {
                zoomView: true,
                minZoom: minScale,
                maxZoom: maxScale
            }
        });

        // Remove duplicate listeners.
        window.network.off("zoom");
        window.network.off("dragEnd");
    };

    const resetView = () => {
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
    };

    const checkScaleAndPosition = () => {
        const currentScale = window.network.getScale();
        const currentPosition = window.network.getViewPosition();

        // If the user zooms out too far, recenter the view.
        if (currentScale < minScale && !isAnimating) {
            resetView();
            return;
        }
        
        // Only check for dragging away from center if we're at or very near the default zoom.
        if (Math.abs(currentScale - minScale) < 0.01 && !isAnimating) {
            const deltaX = Math.abs(currentPosition.x - initialPosition.x);
            const deltaY = Math.abs(currentPosition.y - initialPosition.y);
            if (deltaX > positionThreshold || deltaY > positionThreshold) {
                resetView();
            }
        }
        // When zoomed in (currentScale > minScale), allow the user to pan without automatic recentering.
    };

    const initInterval = setInterval(() => {
        if (window.network && typeof window.network.getScale === "function") {
            setupZoom();
            // Listen to zoom and drag events to monitor the scale and position.
            window.network.on("zoom", checkScaleAndPosition);
            window.network.on("dragEnd", checkScaleAndPosition);
            // Also update zoom bounds on window resize.
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
        "borderWidth": 2,
        "shape": "dot",
        "size":20,
        "font": {
            "size": 12,
            "color": "#000000",
            "face": "arial"
        }
    },
    "edges": {
        "width": 1,
        "color": {
            "opacity": 0.6
        },
        "smooth": {
            "type": "continuous",
            "forceDirection": "none"
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

net.save_graph("graph.html")
HtmlFile = open("graph.html", "r", encoding="utf-8")
source_code = HtmlFile.read()

st.components.v1.html(
    source_code + toggle_script + zoom_control_script,
    height=1000,
    scrolling=True,
)
