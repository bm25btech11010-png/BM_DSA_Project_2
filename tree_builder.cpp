
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cctype>
using namespace std;

// Dynamic Programming

int editDistance(const string &a, const string &b)
{
    int m = a.size(), n = b.size();

    vector<vector<int>> dp(m + 1, vector<int>(n + 1));

    for (int i = 0; i <= m; i++)
        dp[i][0] = i;

    for (int j = 0; j <= n; j++)
        dp[0][j] = j;

    for (int i = 1; i <= m; i++)
    {
        for (int j = 1; j <= n; j++)
        {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;

            dp[i][j] = min({
                dp[i - 1][j] + 1,
                dp[i][j - 1] + 1,
                dp[i - 1][j - 1] + cost
            });
        }
    }

    return dp[m][n];
}

// Graph Theory

struct Edge
{
    int u, v, weight;
};

struct DSU
{
    vector<int> parent, rank;

    DSU(int n)
    {
        parent.resize(n);
        rank.assign(n, 0);

        for (int i = 0; i < n; i++)
            parent[i] = i;
    }

    int find(int x)
    {
        if (parent[x] != x)
            parent[x] = find(parent[x]);

        return parent[x];
    }

    bool unite(int a, int b)
    {
        a = find(a);
        b = find(b);

        if (a == b)
            return false;

        if (rank[a] < rank[b])
            swap(a, b);

        parent[b] = a;

        if (rank[a] == rank[b])
            rank[a]++;

        return true;
    }
};

void kruskal(int n, const vector<vector<int>> &dist,
             const vector<string> &names)
{
    vector<Edge> edges;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            edges.push_back({i, j, dist[i][j]});
        }
    }

    sort(edges.begin(), edges.end(),
         [](const Edge &a, const Edge &b)
         {
             return a.weight < b.weight;
         });

    DSU dsu(n);
    int totalWeight = 0;

    cout << "\n========== GRAPH THEORY ==========\n";
    cout << "Kruskal's Minimum Spanning Tree:\n\n";

    for (const auto &e : edges)
    {
        if (dsu.unite(e.u, e.v))
        {
            cout << names[e.u] << " -- "
                 << names[e.v] << " : "
                 << e.weight << "\n";

            totalWeight += e.weight;
        }
    }

    cout << "\nTotal MST Weight: " << totalWeight << "\n";
}

// Greedy Algorithm

struct TreeNode
{
    string name;
    TreeNode *left = nullptr;
    TreeNode *right = nullptr;
    double height = 0;
};

struct Cluster
{
    vector<int> members;
    TreeNode *node;
};

string newick(TreeNode *node)
{
    if (!node->left && !node->right)
        return node->name;

    double leftBranch = node->height - node->left->height;
    double rightBranch = node->height - node->right->height;

    return "(" + newick(node->left) + ":" +
           to_string(max(0.0, leftBranch)) + "," +
           newick(node->right) + ":" +
           to_string(max(0.0, rightBranch)) + ")";
}

void upgma(int n, const vector<vector<int>> &dist,
           const vector<string> &names)
{
    vector<Cluster> clusters;

    for (int i = 0; i < n; i++)
    {
        TreeNode *leaf = new TreeNode();
        leaf->name = names[i];

        clusters.push_back({{i}, leaf});
    }

    cout << "\n========== GREEDY ALGORITHM ==========\n";
    cout << "UPGMA Cluster Merging Steps:\n\n";

    while (clusters.size() > 1)
    {
        double best = 1e18;
        int a = -1, b = -1;

        for (int i = 0; i < (int)clusters.size(); i++)
        {
            for (int j = i + 1; j < (int)clusters.size(); j++)
            {
                double sum = 0;
                int count = 0;

                for (int x : clusters[i].members)
                {
                    for (int y : clusters[j].members)
                    {
                        sum += dist[x][y];
                        count++;
                    }
                }

                double avg = sum / count;

                if (avg < best)
                {
                    best = avg;
                    a = i;
                    b = j;
                }
            }
        }

        cout << "Merging clusters at distance: "
             << fixed << setprecision(2) << best << "\n";

        TreeNode *parent = new TreeNode();
        parent->left = clusters[a].node;
        parent->right = clusters[b].node;
        parent->height = best / 2.0;

        Cluster merged;
        merged.node = parent;
        merged.members = clusters[a].members;

        for (int x : clusters[b].members)
            merged.members.push_back(x);

        clusters.erase(clusters.begin() + b);
        clusters.erase(clusters.begin() + a);
        clusters.push_back(merged);
    }

    cout << "\nUPGMA Phylogenetic Tree (Newick Format):\n";
    cout << fixed << setprecision(2);
    cout << newick(clusters[0].node) << ";\n";
}

// ==========================================
// MAIN FUNCTION
// ==========================================

int main()
{
    int n;

    cout << "Enter number of DNA sequences (>= 2): ";
    cin >> n;

    vector<string> names(n), sequences(n);

    cout << "\nEnter organism names and DNA sequences "
         << "(use A, T, G, C):\n";

    for (int i = 0; i < n; i++)
    {
        cin >> names[i] >> sequences[i];

        for (char &c : sequences[i])
        {
            c = toupper(static_cast<unsigned char>(c));
        }
    }

    vector<vector<int>> dist(n, vector<int>(n, 0));

    cout << "\n========== DYNAMIC PROGRAMMING ==========\n";
    cout << "Pairwise DNA Edit Distance Matrix:\n\n";

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            dist[i][j] = dist[j][i] =
                editDistance(sequences[i], sequences[j]);
        }
    }

    cout << setw(12) << " ";

    for (const auto &name : names)
        cout << setw(12) << name;

    cout << "\n";

    for (int i = 0; i < n; i++)
    {
        cout << setw(12) << names[i];

        for (int j = 0; j < n; j++)
            cout << setw(12) << dist[i][j];

        cout << "\n";
    }

    kruskal(n, dist, names);
    
    upgma(n, dist, names);

    return 0;
}
