#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

/**
 * f(x, y, z) = x*y + z(x+y)
 * df = A dx + B dy + C dz
 */

enum ArithmeticOp { Add, Minus, Multiply, Divide, Leaf };

struct TreeNode {
    TreeNode(const string &label_, ArithmeticOp op_, double value_) : label{label_}, op{op_}, value{value_} {
        set_variable_name();
        value_evaluated = true;
    }
    TreeNode(const string &label_, ArithmeticOp op_) : label{label_}, op{op_} { set_variable_name(); }

    void set_variable_name() {
        if (label.find(' ') != string::npos) {
            variable_name = label.substr(0, label.find(' '));
        } else {
            variable_name = label;
        }
    }

    TreeNode *left{nullptr}, *right{nullptr};
    string variable_name;
    string label;
    ArithmeticOp op{Leaf};
    double value{0};
    double partial_derivative{0};
    bool value_evaluated{false};
};

inline ostream &operator<<(ostream &s, const TreeNode &data) {
    s << data.variable_name << " = " << data.value << "\t | partial derivative = " << data.partial_derivative << "\t | "
      << data.label;
    return s;
}

class AutoDiffTree {
public:
    AutoDiffTree() {}

    void set_root(const string &label_, ArithmeticOp op) {
        root = new TreeNode(label_, op);
        root->partial_derivative = 1;
    }

    void in_order() {  // 中序遍历
        in_order_recursive(root);
    }

    void pre_order() {  // 前序遍历
        pre_order_recursive(root);
    };

    void post_order() {  // 后续遍历
        post_order_recursive(root);
    };

    void print(TreeNode *currentNode);
    void forward_evaluation(TreeNode *currentNode);
    void reverse_gradient(TreeNode *currentNode);

private:
    void in_order_recursive(TreeNode *currentNode);
    void pre_order_recursive(TreeNode *currentNode);
    void post_order_recursive(TreeNode *currentNode);

public:
    TreeNode *root{nullptr};
};

void AutoDiffTree::in_order_recursive(TreeNode *currentNode) {
    if (currentNode == nullptr) return;
    in_order_recursive(currentNode->left);

    in_order_recursive(currentNode->right);
}

void AutoDiffTree::pre_order_recursive(TreeNode *currentNode) {
    if (currentNode == nullptr) return;

    pre_order_recursive(currentNode->left);
    pre_order_recursive(currentNode->right);
}

void AutoDiffTree::post_order_recursive(TreeNode *currentNode) {
    if (currentNode == nullptr) return;
    post_order_recursive(currentNode->left);
    post_order_recursive(currentNode->right);
}

void AutoDiffTree::print(TreeNode *currentNode) {
    if (currentNode == nullptr) return;
    print(currentNode->left);
    print(currentNode->right);
    cout << *currentNode << endl;
}

void AutoDiffTree::forward_evaluation(TreeNode *currentNode) {
    if (currentNode == nullptr || currentNode->value_evaluated) return;
    forward_evaluation(currentNode->left);
    forward_evaluation(currentNode->right);

    switch (currentNode->op) {
        case Add: {
            currentNode->value = currentNode->left->value + currentNode->right->value;
            break;
        }
        case Minus: {
            currentNode->value = currentNode->left->value - currentNode->right->value;
            break;
        }
        case Multiply: {
            currentNode->value = currentNode->left->value * currentNode->right->value;
            break;
        }
        case Divide: {
            currentNode->value = currentNode->left->value / currentNode->right->value;
            break;
        }
        default:
            break;
    }
}

void AutoDiffTree::reverse_gradient(TreeNode *currentNode) {
    if (currentNode == nullptr || currentNode->left == nullptr || currentNode->right == nullptr) return;

    switch (currentNode->op) {
        case Add: {
            currentNode->left->partial_derivative += 1 * currentNode->partial_derivative;
            currentNode->right->partial_derivative += 1 * currentNode->partial_derivative;
            break;
        }
        case Minus: {
            currentNode->left->partial_derivative += 1 * currentNode->partial_derivative;
            currentNode->right->partial_derivative += -1 * currentNode->partial_derivative;
            break;
        }
        case Multiply: {
            currentNode->left->partial_derivative += currentNode->right->value * currentNode->partial_derivative;
            currentNode->right->partial_derivative += currentNode->left->value * currentNode->partial_derivative;
            break;
        }
        case Divide: {
            currentNode->left->partial_derivative +=
                (1.0 / currentNode->right->value) * currentNode->partial_derivative;
            currentNode->right->partial_derivative +=
                (-currentNode->left->value / std::pow(currentNode->right->value, 2)) * currentNode->partial_derivative;
            break;
        }
        default:
            break;
    }
    reverse_gradient(currentNode->left);
    reverse_gradient(currentNode->right);
}

int main() {
    AutoDiffTree tree;
    tree.set_root("w4 = w1 + w3", Add);
    tree.root->right = new TreeNode("w1 = x * y", Multiply);
    tree.root->left = new TreeNode("w3 = w2 * z", Multiply);
    tree.root->left->left = new TreeNode("w2 = x * y", Add);
    tree.root->left->right = new TreeNode("z", Leaf, 1);
    auto *x = new TreeNode("x", Leaf, 1);
    auto *y = new TreeNode("y", Leaf, 1);
    tree.root->left->left->left = x;
    tree.root->left->left->right = y;
    tree.root->right->left = x;
    tree.root->right->right = y;

    cout << "init tree:\n";
    tree.print(tree.root);

    tree.forward_evaluation(tree.root);
    cout << "\nafter forward evaluation:\n";
    tree.print(tree.root);

    tree.reverse_gradient(tree.root);
    cout << "\nafter reverse gradient:\n";
    tree.print(tree.root);

    cout << "\nend\n\n";
    return 0;
}
