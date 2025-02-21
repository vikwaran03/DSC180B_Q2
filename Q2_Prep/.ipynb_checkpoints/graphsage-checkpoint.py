import torch
import torch_geometric
from torch_geometric import utils
from torch_geometric.nn import SAGEConv
import torch.nn as nn
import torch.nn.functional as F

from scipy.sparse import coo_matrix
import pandas as pd
import numpy as np

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Load in Hi-c Matrices and node information, build graph
ec_hic = np.load('data/ec_adj_mat_t25.npy')
hsr_hic = np.load('data/hsr_adj_mat_t25.npy')

ec_df = pd.read_csv('data/ec_cleaned.csv')
hsr_df = pd.read_csv('data/hsr_cleaned.csv')

hsr_feats = torch.tensor(hsr_df[['read_count', 'total_genes']].to_numpy())
hsr_labels = torch.zeros(hsr_feats.shape[0])

ec_feats = torch.tensor(ec_df[['read_count', 'total_genes']].to_numpy())
ec_labels = torch.ones(ec_feats.shape[0])

def hic_to_sparse(hic_mat):
    adj_mat = np.triu(hic_mat)
    sparse_adj = coo_matrix(adj_mat)

    return utils.from_scipy_sparse_matrix(sparse_adj)

hsr_edge_index, hsr_edge_attr = hic_to_sparse(hsr_hic)
hsr_graph = torch_geometric.data.Data(edge_index = hsr_edge_index, edge_attr = hsr_edge_attr, x = hsr_feats, y = hsr_labels)

ec_edge_index, ec_edge_attr = hic_to_sparse(ec_hic)
ec_graph = torch_geometric.data.Data(edge_index = ec_edge_index, edge_attr = ec_edge_attr, x = ec_feats, y = ec_labels)

x = torch.cat([ec_feats, hsr_feats], dim=0)
hsr_edge_index = hsr_edge_index + ec_labels.shape[0]
edge_index = torch.cat([ec_edge_index, hsr_edge_index], dim=1)
edge_attr = torch.cat([ec_edge_attr, hsr_edge_attr], dim=0)
labels = torch.cat([ec_labels, hsr_labels], dim=0)

G = torch_geometric.data.Data(edge_index = edge_index, edge_attr = edge_attr, x = x, y = labels)

# Define GraphSAGE
class GraphSAGE(nn.Module):
    def __init__(self, num_feat, graph_conv_layer_sizes, lin_hidden_sizes, num_classes):
        super().__init__()
        self.num_feat = num_feat
        self.num_classes = num_classes
        self.embeddings = None

        self.conv1 = SAGEConv(self.num_feat, graph_conv_layer_sizes[0])
        self.conv2 = SAGEConv(graph_conv_layer_sizes[0], graph_conv_layer_sizes[1])
        
        self.lin1 = nn.Linear(lin_hidden_sizes[0], lin_hidden_sizes[1])
        self.lin2 = nn.Linear(lin_hidden_sizes[1], self.num_classes)
        
        self.loss_calc = nn.CrossEntropyLoss()
        self.torch_softmax = nn.Softmax(dim=1)

    def forward(self, data):
        x, edge_index, edge_attr = data.x, data.edge_index, data.edge_attr

        ### Graph convolution module
        h = self.conv1(x, edge_index)
        h = torch.relu(h)
        h = self.conv2(h, edge_index)
        h = torch.relu(h)
            
        self.embeddings = h
        scores = h
        
        ### Linear module
        scores = self.lin1(scores)
        scores = torch.relu(scores)
        scores = self.lin2(scores)
        
        return scores

    def loss(self, scores, labels):
        xent_loss = self.loss_calc(scores, labels)
        return xent_loss
    
    def calc_softmax_pred(self, scores):
        softmax = self.torch_softmax(scores)
        predicted = torch.argmax(softmax, 1)
        return softmax, predicted

    def get_embeddings(self):
        if self.embeddings is not None:
            return self.embeddings

        print("Untrained Model: Please train model first")
        return None

# separate train/test sets
def split_data(G, train_ratio=0.7, val_ratio=0.15, test_ratio=0.15, seed=42):
    assert train_ratio + val_ratio + test_ratio == 1, "Ratios must sum to 1."
    
    torch.manual_seed(seed)
    num_nodes = G.x.shape[0]  # Total number of nodes
    indices = torch.randperm(num_nodes)  # Shuffle node indices randomly
    
    train_size = int(train_ratio * num_nodes)
    val_size = int(val_ratio * num_nodes)
    train_idx = indices[:train_size]
    val_idx = indices[train_size:train_size + val_size]
    test_idx = indices[train_size + val_size:]

    # Create boolean masks
    train_mask = torch.zeros(num_nodes, dtype=torch.bool)
    val_mask = torch.zeros(num_nodes, dtype=torch.bool)
    test_mask = torch.zeros(num_nodes, dtype=torch.bool)

    train_mask[train_idx] = True
    val_mask[val_idx] = True
    test_mask[test_idx] = True

    # Attach masks to the graph
    G.train_mask = train_mask
    G.val_mask = val_mask
    G.test_mask = test_mask

    return G

# Train and Test Functions
def train(model, optimizer, graph, device):
    model.train()
    optimizer.zero_grad()
    scores = model(graph)
    loss = model.loss(scores[graph.train_mask], graph.y[graph.train_mask])
    loss.backward()
    optimizer.step()
    
    return loss.item(), scores

def test(model, graph, device):
    model.eval()
    scores = model(graph)
    softmax, predicted = model.calc_softmax_pred(scores)
    accs = []
    losses = []
    
    for mask in [graph.train_mask, graph.val_mask, graph.test_mask]:
        loss = model.loss(scores[mask], graph.y[mask]).item()
        correct = (predicted[mask] == graph.y[mask]).sum().item()
        acc = correct / mask.sum().item()
        accs.append(acc)
        losses.append(loss)
    
    return accs, losses

# Model inputs + layer information
num_features = G.num_node_features
num_classes = 2
graph_sage_layer_sizes = [8,16]
linear_layer_sizes = [16,8]

# Hyperparameters
learning_rate = 0.01
weight_decay = 5e-4
num_epochs = 5000
patience = 50  # Number of epochs with no improvement before stopping
min_delta = 0.001  # Minimum change to qualify as an improvement

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
best_models = []
all_metrics = []

# Training loop - Entire graph in one model
for _ in range(5):
    # Initialize the model, optimizer, and send to GPU
    model = GraphSAGE(num_features, graph_sage_layer_sizes, linear_layer_sizes, num_classes).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)

    # Assign train and test masks to graph
    G.x = G.x.float()
    G.y = G.y.long()
    G = split_data(G)
    G = G.to(device)
    
    train_accs = []
    val_accs = []
    test_accs = []
    train_losses = []
    val_losses = []
    test_losses = []
    
    best_val_loss = float('inf')
    patience_counter = 0
    
    best_model = None
    
    for epoch in range(1, num_epochs + 1):
        loss, score = train(model, optimizer, G, device)
        (train_acc, val_acc, test_acc), (train_loss, val_loss, test_loss) = test(model, G, device)
    
        # Track metrics
        train_losses.append(loss)
        valid_losses.append(val_loss)
        test_losses.append(test_loss)
        train_accs.append(train_acc)
        val_accs.append(val_acc)
        test_accs.append(test_acc)
    
        if epoch % 500 == 0:
            print(f"Epoch {epoch:03d}, Loss: {loss:.4f}, Train Acc: {train_acc:.4f}, Val Acc: {val_acc:.4f}, Test Acc: {test_acc:.4f}")
    
        # Early stopping check - only start checking after 3000 epochs
        if epoch >= 3000:
            if val_loss < best_val_loss - min_delta:
                best_val_loss = val_loss
                patience_counter = 0
                best_model = model.state_dict()  # Save best model
            else:
                patience_counter += 1

            if epoch == 5000 and patience_counter < patience:
                best_models.append(best_model)
                all_metrics.append(np.array([train_accs, test_accs, val_accs, train_losses, valid_losses, test_losses]))
                
            # Stop if validation loss doesn't improve for 'patience' epochs
            if patience_counter >= patience:
                print(f"Early stopping triggered at epoch {epoch}. Best validation loss: {best_val_loss:.4f}")
                best_models.append(best_model)
                all_metrics.append(np.array([train_accs, test_accs, val_accs, train_losses, valid_losses, test_losses]))
                break

# save model - figure out how to save best performing model from the 5 runs
#torch.save(save_model.state_dict(), "graphsage_adj.pth")

# Plot PCA of embeddings
embeddings = model.get_embeddings().cpu().detach().numpy()
pca = PCA(n_components=2)
pca_embeds = pca.fit_transform(embeddings)
ec_mask = (G.y == 1).cpu()

plt.figure(figsize=(8,8))
plt.scatter(pca_embeds[:, 0][ec_mask], pca_embeds[:, 1][ec_mask], color='red', label='ecDNA', alpha=0.7)
plt.scatter(pca_embeds[:, 0][~ec_mask], pca_embeds[:, 1][~ec_mask], color='blue', label='HSR', alpha=0.7)
plt.title("GraphSAGE Embeddings PCA")
plt.ylim(-110, 110)
plt.legend()
plt.savefig('sage_embeddings.png')

# Plot test accuracy and loss over epochs
plt.figure(figsize=(8, 8))
plt.plot(test_accs)
plt.xlabel("Epochs (x50)")
plt.ylabel("Test Accuracy")
plt.title("Test Accuracy Across GraphSAGE Run")
plt.legend()
plt.savefig('sage_acc.png')

plt.figure(figsize=(8, 8))
plt.plot(losses)
plt.xlabel("Epochs (x50)")
plt.ylabel("Loss")
plt.title("Loss Across GraphSAGE Run")
plt.legend()
plt.savefig('sage_loss.png')