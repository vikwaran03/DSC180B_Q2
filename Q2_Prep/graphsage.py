import torch
import torch_geometric
from torch_geometric import utils
from torch_geometric.nn import SAGEConv
import torch.nn as nn
import torch.nn.functional as F
from torch.optim.lr_scheduler import ReduceLROnPlateau

from scipy.sparse import coo_matrix
import pandas as pd
import numpy as np
import os

from sklearn.metrics import precision_score, recall_score, f1_score, confusion_matrix
from sklearn.decomposition import PCA 
import seaborn as sns
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
        self.embeddings = None
        self.conv1 = SAGEConv(num_feat, graph_conv_layer_sizes[0])
        self.conv2 = SAGEConv(graph_conv_layer_sizes[0], graph_conv_layer_sizes[1])

        self.lin1 = nn.Linear(lin_hidden_sizes[0], lin_hidden_sizes[1])
        self.lin2 = nn.Linear(lin_hidden_sizes[1], num_classes)
            
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
    precisions = []
    recalls = []
    f1_scores = []

    for mask in [graph.train_mask, graph.val_mask, graph.test_mask]:
        loss = model.loss(scores[mask], graph.y[mask]).item()
        correct = (predicted[mask] == graph.y[mask]).sum().item()
        acc = correct / mask.sum().item()
        accs.append(acc)
        losses.append(loss)

        # Convert tensors to numpy arrays for sklearn
        y_true = graph.y[mask].cpu().numpy()
        y_pred = predicted[mask].cpu().numpy()

        # Want high recall. High recall = fewer FN. FN is worse than FP. Predicting no cancer when there is cancer is bad
        precision = precision_score(y_true, y_pred, average='macro', zero_division=0)
        recall = recall_score(y_true, y_pred, average='macro', zero_division=0)
        f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)

        precisions.append(precision)
        recalls.append(recall)
        f1_scores.append(f1)

    return accs, losses, precisions, recalls, f1_scores, y_true, y_pred

# Model inputs + layer information
num_features = G.num_node_features
num_classes = 2
graph_sage_layer_sizes = [8,16]
linear_layer_sizes = [16,8]

# Hyperparameters
learning_rate = 0.01
weight_decay = 5e-4
num_epochs = 5000
patience = 100  # Number of epochs with no improvement before stopping
min_delta = 0.001  # Minimum change to qualify as an improvement

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
# Assign train and test masks to graph
G.x = G.x.float()
G.y = G.y.long()
G = split_data(G)
G = G.to(device)

best_models = []
all_metrics = []
all_preds = []

# Training loop - Entire graph in one model
for _ in range(5):
    # Initialize the model, optimizer, and send to GPU
    model = GraphSAGE(num_features, graph_sage_layer_sizes, linear_layer_sizes, num_classes).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    scheduler = ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=patience//2, verbose=True) #Reduces LR by 50% when training plateaus
    
    train_accs = []
    val_accs = []
    test_accs = []
    train_losses = []
    val_losses = []
    test_losses = []
    train_precisions = []
    val_precisions = []
    test_precisions = []
    train_recalls = []
    val_recalls = []
    test_recalls = []
    train_f1_scores = []
    val_f1_scores = []
    test_f1_scores = []
    preds = []
    
    best_val_loss = float('inf')
    patience_counter = 0
    
    best_model = None
    
    for epoch in range(1, num_epochs + 1):
        loss, score = train(model, optimizer, G, device)
        accs, losses, precisions, recalls, f1_scores, y_true, y_pred = test(model, G, device)
        #(train_acc, val_acc, test_acc), (train_loss, val_loss, test_loss) = test(model, G, device)
        train_acc = accs[0]
        val_acc = accs[1]
        test_acc = accs[2]
        val_loss = losses[1]
    
        train_losses.append(loss)
        val_losses.append(val_loss)
        test_losses.append(losses[2])
        train_accs.append(train_acc)
        val_accs.append(val_acc)
        test_accs.append(test_acc)
        train_precisions.append(precisions[0])
        val_precisions.append(precisions[1])
        test_precisions.append(precisions[2])
        train_recalls.append(recalls[0])
        val_recalls.append(recalls[1])
        test_recalls.append(recalls[2])
        train_f1_scores.append(f1_scores[0])
        val_f1_scores.append(f1_scores[1])
        test_f1_scores.append(f1_scores[2])
        preds.append([y_true, y_pred])
    
        if epoch % 500 == 0:
            print(f"Epoch {epoch:03d}, Loss: {loss:.4f}, Train Acc: {train_acc:.4f}, Val Acc: {val_acc:.4f}, Test Acc: {test_acc:.4f}")
    
        # Early stopping check - only start checking after 3500 epochs
        if epoch >= 3500:
            if val_loss < best_val_loss - min_delta:
                best_val_loss = val_loss
                patience_counter = 0
                best_model = model  # Save best model
            else:
                patience_counter += 1

            # Adjust learning rate
            scheduler.step(val_loss)

            if epoch == 5000 and patience_counter < patience:
                best_models.append(best_model)
                all_metrics.append(np.array([train_accs, test_accs, val_accs, train_losses, val_losses, test_losses]))
                all_preds.append(preds)
                
            # Stop if validation loss doesn't improve for 'patience' epochs
            if patience_counter >= patience:
                print(f"Early stopping triggered at epoch {epoch}. Best validation loss: {best_val_loss:.4f}")
                best_models.append(best_model)
                all_metrics.append(np.array([train_accs, test_accs, val_accs, train_losses, val_losses, test_losses, train_precisions, val_precisions, test_precisions, train_recalls, val_recalls, test_recalls, train_f1_scores, val_f1_scores, test_f1_scores]))
                all_preds.append(preds)
                break

# save final results from all runs and find best performing model
train_accs_across_runs = [all_metrics[i][0][-patience] for i in range(len(all_metrics))]
val_accs_across_runs = [all_metrics[i][2][-patience] for i in range(len(all_metrics))]
test_accs_across_runs = [all_metrics[i][1][-patience] for i in range(len(all_metrics))]
test_precision_across_runs = [all_metrics[i][8][-patience] for i in range(len(all_metrics))]
test_recall_across_runs = [all_metrics[i][11][-patience] for i in range(len(all_metrics))]
test_f1score_across_runs = [all_metrics[i][14][-patience] for i in range(len(all_metrics))]
max_acc = 0
best_model_idx = None
for i,acc in enumerate(test_accs_across_runs):
    if acc > max_acc:
        max_acc = acc
        best_model_idx = i

if not os.path.exists('images'):
    os.makedirs('images')

if not os.path.exists('model'):
    os.makedirs('model')

if not os.path.exists('outputs'):
    os.makedirs('outputs')

res_table_runs = pd.DataFrame({'Train Accuracy': train_accs_across_runs, 'Valdiation Accuracy': val_accs_across_runs, 'Test Accuracy': test_accs_across_runs, 'Test Precision': test_precision_across_runs, 'Test Recall': test_recall_across_runs, 'Test F1 Score': test_f1score_across_runs}, index = ['Run 1', 'Run 2', 'Run 3', 'Run 4', 'Run 5'])
res_table_runs.to_csv('outputs/sage_run_results_table.csv')

torch.save(best_models[best_model_idx].state_dict(), "model/graphsage_adj.pth")
best_metrics = all_metrics[best_model_idx]
best_preds = all_preds[best_model_idx][-patience]

best_model_table = pd.DataFrame({'Train': [best_metrics[0][-patience], best_metrics[6][-patience], best_metrics[9][-patience], best_metrics[12][-patience], best_metrics[3][-patience]], 'Validation': [best_metrics[1][-patience], best_metrics[7][-patience], best_metrics[10][-patience], best_metrics[13][-patience], best_metrics[4][-patience]], 'Test': [best_metrics[2][-patience], best_metrics[8][-patience], best_metrics[11][-patience], best_metrics[14][-patience], best_metrics[5][-patience]]}, index = ['Accuracy', 'Precision', 'Recall', 'F1 Score', 'Loss'])
best_model_table.to_csv('outputs/sage_best_model_table.csv')

# Plot test accuracy and loss across runs
plt.figure(figsize=(8, 8))
for i, row in enumerate([all_metrics[i][1] for i in range(len(all_metrics))]):
    plt.plot(row, label=f"Run {i+1}")
plt.xlabel("Epochs")
plt.ylabel("Test Accuracy")
plt.title("Test Accuracy Across 5 GraphSAGE Runs")
plt.legend()
plt.savefig('images/sage_acc.png')

plt.figure(figsize=(8, 8))
for i, row in enumerate([all_metrics[i][3] for i in range(len(all_metrics))]):
    plt.plot(row, label=f"Run {i+1}")
plt.xlabel("Epochs")
plt.ylabel("Loss")
plt.ylim(0,10)
plt.title("Training Loss Across 5 GraphSAGE Runs")
plt.legend()
plt.savefig('images/sage_loss.png')

# Plot PCA, Accuracies, Loss for best model
# Plot PCA of embeddings
embeddings = best_models[best_model_idx].get_embeddings().cpu().detach().numpy()
np.save("outputs/sage_embeddings.npy", embeddings)
pca = PCA(n_components=2)
pca_embeds = pca.fit_transform(embeddings)
ec_mask = (G.y == 1).cpu()
np.save("outputs/ec_pca.npy", pca_embeds[ec_mask])
np.save("outputs/hsr_pca.npy", pca_embeds[~ec_mask])

plt.figure(figsize=(8,8))
plt.scatter(pca_embeds[:, 0][ec_mask], pca_embeds[:, 1][ec_mask], color='red', label='ecDNA', alpha=0.7)
plt.scatter(pca_embeds[:, 0][~ec_mask], pca_embeds[:, 1][~ec_mask], color='blue', label='HSR', alpha=0.7)
plt.title("GraphSAGE Embeddings PCA")
plt.ylim(-110, 110)
plt.legend()
plt.savefig('images/sage_embeddings.png')

plt.figure(figsize=(8, 8))
plt.plot(best_metrics[0], label="Train Accuracy")
plt.plot(best_metrics[1], label="Test Accuracy")
plt.plot(best_metrics[2], label="Validation Accuracy")
plt.xlabel("Epochs")
plt.ylabel("Accuracy")
plt.title("Accuracy Curves Best Model")
plt.legend()
plt.savefig('images/sage_best_model_acc_curves.png')

plt.figure(figsize=(8, 8))
plt.plot(best_metrics[3], label="Train Loss")
plt.plot(best_metrics[4], label="Validation Loss")
plt.plot(best_metrics[5], label="Test Loss")
plt.xlabel("Epochs")
plt.ylabel("Loss")
plt.ylim(0,10)
plt.title("Loss Curves Best Model")
plt.legend()
plt.savefig('images/sage_best_model_loss_curves.png')

cm = confusion_matrix(best_preds[0], best_preds[1])
plt.figure(figsize=(8, 8))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['HSR', 'ecDNA'], yticklabels=['HSR', 'ecDNA'])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix on Test Set')
plt.savefig('images/sage_confusion_matrix.png')