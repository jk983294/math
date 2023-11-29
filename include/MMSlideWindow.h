#ifndef ORNATE_MMSLIDEWINDOW_H
#define ORNATE_MMSLIDEWINDOW_H

#include <iostream>
#include <vector>

/**
 * a sliding window size K, keep track of min max value within window
 */

enum MmswColor { BLACK, RED };

template <typename T>
struct MmswRBNode {
    T key;        // data
    int index;    // index of the data in the pool
    size_t size;  // size of the sub tree
    MmswColor color;
    MmswRBNode* left;
    MmswRBNode* right;
    MmswRBNode* parent;
};

template <typename T>
struct MaxMinSlideWindow {
    explicit MaxMinSlideWindow(size_t length_k_) : length_k{length_k_} {
        create_nil();
        nodes.resize(length_k_);
        pool = new MmswRBNode<T>[length_k_ + 1];
        pool_root = pool;
        for (size_t i = 0; i < length_k_; i++) {
            pool[i].right = &pool[i + 1];
            pool[i + 1].parent = &pool[i];  // gist of length_k + 1 pool size
        }
    }

    // copy constructor
    MaxMinSlideWindow(const MaxMinSlideWindow& mmw) : nil_key(mmw.nil_key), index(mmw.index), length_k(mmw.length_k) {
        create_nil();
        nodes.resize(length_k);
        pool = new MmswRBNode<T>[length_k + 1];
        pool_root = pool;
        for (int i = 0; i < length_k; i++) {
            pool[i].right = &pool[i + 1];
            pool[i + 1].parent = &pool[i];
        }
        if (mmw.root != mmw.nil) {
            root = pool_pop();
            root->parent = nil;
            copy_tree(root, mmw.root, mmw);
        } else {
            root = nil;
        }
    }

    void create_nil() {
        if (nil == nullptr) nil = new MmswRBNode<T>;
        nil->color = BLACK;
        nil->size = 0;
        nil->key = this->nil_key;
        root = nil;
    }

    MaxMinSlideWindow& operator=(const MaxMinSlideWindow& mmw) {
        nil_key = mmw.nil_key;
        index = mmw.index;
        length_k = mmw.length_k;
        create_nil();
        nodes.resize(length_k);
        pool = new MmswRBNode<T>[length_k + 1];
        pool_root = pool;
        for (int i = 0; i < length_k; i++) {
            pool[i].right = &pool[i + 1];
            pool[i + 1].parent = &pool[i];
        }
        if (mmw.root != mmw.nil) {
            root = pool_pop();
            root->parent = nil;
            copy_tree(root, mmw.root, mmw);
        } else {
            root = nil;
        }
        return (*this);
    }

    ~MaxMinSlideWindow();

    void clear();

    /**
     * define nan(not a number)
     */
    void set_nan(T nan) { nil_key = nan; }

    /**
     *  main operation of MaxMinSlideWindow
     *  1. insert a data into data window.
     *  2. remove oldest data if data window is full.
     *  3. maintain the tree depth of the red-black-tree.
     */
    void add_new(const T& x);

    // get constant reference of kth min value in data window, k >= 0
    inline const T& mink_val(size_t k) const {
        k++;
        if (k > root->size) {
            return (nil_key);
        }
        return (Select(root, k)->key);
    }

    // get constant reference of kth max value in data window, k >= 0
    inline const T& maxk_val(size_t k) const {
        k++;
        if (k > root->size) {
            return (nil_key);
        }
        if (index < length_k)
            return (Select(root, index - k + 1)->key);
        else
            return (Select(root, length_k - k + 1)->key);
    }

    inline int mink_pos(size_t k) const {
        k++;
        if (k > root->size) {
            return -1;
        } else {
            if (index < length_k) {
                return (Select(root, k)->index);
            } else {
                return static_cast<int>(Select(root, k)->index + length_k - index);
            }
        }
    }

    inline int maxk_pos(size_t k) const {
        k++;
        if (k > root->size) {
            return -1;
        }
        if (index < length_k) {
            return (Select(root, index - k + 1)->index);
        } else {
            return static_cast<int>(Select(root, length_k - k + 1)->index + length_k - index);
        }
    }

    // get constant reference of min value in data window
    inline const T& min_val() const { return (Minimum(root)->key); }

    // get constant reference of max value in data window
    inline const T& max_val() const { return (Maximum(root)->key); }

    inline int min_pos() const {
        if (index < length_k) {
            return (Minimum(root)->index);
        } else {
            return static_cast<int>(Minimum(root)->index - (index - length_k));
        }
    }

    inline int max_pos() const {
        if (index < length_k) {
            return (Maximum(root)->index);
        } else {
            return static_cast<int>((Maximum(root)->index) - index + length_k);
        }
    }

    inline const size_t size() { return root->size; }

    void print_all_data() {
        if (index < length_k) {
            for (size_t i = 0; i < index; ++i) {
                std::cout << nodes[i]->key << " ";
            }
            std::cout << "all data" << std::endl;
        } else {
            for (size_t i = 0; i < length_k; ++i) {
                std::cout << nodes[i]->key << " ";
            }
            std::cout << "all data" << std::endl;
        }
    }

private:
    MmswRBNode<T>* Successor(MmswRBNode<T>* x);
    void InsertNode(MmswRBNode<T>* z);
    void copy_tree(MmswRBNode<T>* node, const MmswRBNode<T>* from, const MaxMinSlideWindow<T>& mmw);
    MmswRBNode<T>* DeleteNode(MmswRBNode<T>* z);
    inline MmswRBNode<T>* Maximum(MmswRBNode<T>* x) const;
    inline MmswRBNode<T>* Minimum(MmswRBNode<T>* x) const;
    inline MmswRBNode<T>* Select(MmswRBNode<T>* x, size_t i) const;  // i is ith element in data window, [1, K]
    inline MmswRBNode<T>* pool_pop();
    inline void pool_push(MmswRBNode<T>* z);

    void LeftRotate(MmswRBNode<T>* x);
    void RightRotate(MmswRBNode<T>* y);
    void InsertFixup(MmswRBNode<T>* z);
    void DeleteFixup(MmswRBNode<T>* x);

    std::vector<MmswRBNode<T>*> nodes;  // pointers point the last length_k Nodes
    MmswRBNode<T>* pool_root{nullptr};  // next available node in pool
    MmswRBNode<T>* root{nullptr};       // root of the tree
    MmswRBNode<T>* nil{nullptr};
    MmswRBNode<T>* pool{nullptr};
    T nil_key;
    size_t index{0};  // The count of the flow after last clear()
    size_t length_k;
};

template <typename T>
void MaxMinSlideWindow<T>::copy_tree(MmswRBNode<T>* node, const MmswRBNode<T>* from, const MaxMinSlideWindow<T>& mmw) {
    node->key = from->key;
    node->index = from->index;
    node->size = from->size;
    node->color = from->color;
    nodes[node->index % length_k] = node;  // no collision because we keep last length_k window data

    if (from->left != mmw.nil) {
        MmswRBNode<T>* left_node = pool_pop();
        left_node->parent = node;
        node->left = left_node;
        copy_tree(left_node, from->left, mmw);
    } else {
        node->left = nil;
    }
    if (from->right != mmw.nil) {
        MmswRBNode<T>* right_node = pool_pop();
        right_node->parent = node;
        node->right = right_node;
        copy_tree(right_node, from->right, mmw);
    } else {
        node->right = nil;
    }
}

template <typename T>
MaxMinSlideWindow<T>::~MaxMinSlideWindow<T>() {
    delete nil;
    delete[] pool;
}

template <typename T>
inline MmswRBNode<T>* MaxMinSlideWindow<T>::pool_pop() {
    pool_root = pool_root->right;
    return pool_root->parent;
}

template <typename T>
inline void MaxMinSlideWindow<T>::pool_push(MmswRBNode<T>* z) {
    pool_root->parent = z;
    z->right = pool_root;
    pool_root = z;
}

template <typename T>
void MaxMinSlideWindow<T>::LeftRotate(MmswRBNode<T>* x) {
    MmswRBNode<T>* y = x->right;
    x->right = y->left;
    if (y->left != nil) y->left->parent = x;
    y->parent = x->parent;
    if (x->parent == nil) {
        root = y;
    } else if (x == x->parent->left) {
        x->parent->left = y;
    } else {
        x->parent->right = y;
    }
    y->left = x;
    x->parent = y;
    y->size = x->size;
    x->size = x->left->size + x->right->size + 1;
}

template <typename T>
void MaxMinSlideWindow<T>::RightRotate(MmswRBNode<T>* y) {
    MmswRBNode<T>* x = y->left;
    y->left = x->right;
    if (x->right != nil) x->right->parent = y;
    x->parent = y->parent;
    if (y->parent == nil) {
        root = x;
    } else if (y == y->parent->right) {
        y->parent->right = x;
    } else {
        y->parent->left = x;
    }
    x->right = y;
    y->parent = x;
    x->size = y->size;
    y->size = y->right->size + y->left->size + 1;
}

template <typename T>
void MaxMinSlideWindow<T>::InsertNode(MmswRBNode<T>* z) {
    MmswRBNode<T>*y = nil, *x = root;
    while (x != nil) {
        x->size = x->size + 1;
        y = x;
        if (z->key < x->key) {
            x = x->left;
        } else {
            x = x->right;
        }
    }
    z->parent = y;
    if (y == nil) {
        root = z;
    } else if (z->key < y->key) {
        y->left = z;
    } else {
        y->right = z;
    }
    z->left = nil;
    z->right = nil;
    z->color = RED;
    z->size = 1;
    InsertFixup(z);
}

template <typename T>
void MaxMinSlideWindow<T>::InsertFixup(MmswRBNode<T>* z) {
    MmswRBNode<T>* y;
    while (z->parent->color == RED) {
        if (z->parent == z->parent->parent->left) {
            y = z->parent->parent->right;
            if (y->color == RED) {
                z->parent->color = BLACK;
                y->color = BLACK;
                z->parent->parent->color = RED;
                z = z->parent->parent;
            } else {
                if (z == z->parent->right) {
                    z = z->parent;
                    LeftRotate(z);
                }
                z->parent->color = BLACK;
                z->parent->parent->color = RED;
                RightRotate(z->parent->parent);
            }
        } else {
            y = z->parent->parent->left;
            if (y->color == RED) {
                z->parent->color = BLACK;
                y->color = BLACK;
                z->parent->parent->color = RED;
                z = z->parent->parent;
            } else {
                if (z == z->parent->left) {
                    z = z->parent;
                    RightRotate(z);
                }
                z->parent->color = BLACK;
                z->parent->parent->color = RED;
                LeftRotate(z->parent->parent);
            }
        }
    }
    root->color = BLACK;
}

template <typename T>
MmswRBNode<T>* MaxMinSlideWindow<T>::DeleteNode(MmswRBNode<T>* z) {
    MmswRBNode<T>*y, *x;
    if (z->left == nil || z->right == nil) {
        y = z;
    } else {
        y = Successor(z);
    }
    if (y->left != nil) {
        x = y->left;
    } else {
        x = y->right;
    }
    x->parent = y->parent;

    if (y->parent == nil) {
        root = x;
    } else if (y == y->parent->left) {
        y->parent->left = x;
    } else {
        y->parent->right = x;
    }

    if (y != z) {
        z->key = y->key;
        /// copy y's satellite data to z.
        z->index = y->index;
        nodes[y->index % (length_k)] = z;
    }

    z = y;
    while (z != root) {
        z = z->parent;
        if (z == nil) break;
        z->size = z->size - 1;
    }
    if (y->color == BLACK) {
        DeleteFixup(x);
    }
    return y;
}

template <typename T>
void MaxMinSlideWindow<T>::DeleteFixup(MmswRBNode<T>* x) {
    MmswRBNode<T>* w;

    while (x != root && x->color == BLACK) {
        if (x == x->parent->left) {
            w = x->parent->right;
            if (w->color == RED) {
                w->color = BLACK;
                x->parent->color = RED;
                LeftRotate(x->parent);
                w = x->parent->right;
            }
            if (w->left->color == BLACK && w->right->color == BLACK) {
                w->color = RED;
                x = x->parent;
            } else {
                if (w->right->color == BLACK) {
                    w->left->color = BLACK;
                    w->color = RED;
                    RightRotate(w);
                    w = x->parent->right;
                }
                w->color = x->parent->color;
                x->parent->color = BLACK;
                w->right->color = BLACK;
                LeftRotate(x->parent);
                x = root;
            }
        } else {
            w = x->parent->left;
            if (w->color == RED) {
                w->color = BLACK;
                x->parent->color = RED;
                RightRotate(x->parent);
                w = x->parent->left;
            }
            if (w->right->color == BLACK && w->left->color == BLACK) {
                w->color = RED;
                x = x->parent;
            } else {
                if (w->left->color == BLACK) {
                    w->right->color = BLACK;
                    w->color = RED;
                    LeftRotate(w);
                    w = x->parent->left;
                }

                w->color = x->parent->color;
                x->parent->color = BLACK;
                w->left->color = BLACK;
                RightRotate(x->parent);
                x = root;
            }
        }
    }
    x->color = BLACK;
}

template <typename T>
inline MmswRBNode<T>* MaxMinSlideWindow<T>::Select(MmswRBNode<T>* x, size_t i) const {
    size_t r = x->left->size + 1;  // x's order = left tree size + itself
    if (i > root->size) {
        return nil;
    }
    if (i == r) {
        return x;
    } else if (i < r) {
        return Select(x->left, i);
    } else {
        return Select(x->right, i - r);
    }
}

template <typename T>
MmswRBNode<T>* MaxMinSlideWindow<T>::Successor(MmswRBNode<T>* x) {
    if (x->right != nil) {
        return Minimum(x->right);
    }
    MmswRBNode<T>* y = x->parent;
    while (y != nil && x == y->right) {
        x = y;
        y = y->parent;
    }
    return y;
}

template <typename T>
inline MmswRBNode<T>* MaxMinSlideWindow<T>::Minimum(MmswRBNode<T>* x) const {
    while (x->left != nil) {
        x = x->left;
    }
    return x;
}

template <typename T>
inline MmswRBNode<T>* MaxMinSlideWindow<T>::Maximum(MmswRBNode<T>* x) const {
    while (x->right != nil) {
        x = x->right;
    }
    return x;
}

template <typename T>
void MaxMinSlideWindow<T>::clear() {
    index = 0;
    nil->color = BLACK;
    nil->size = 0;
    root = nil;
    pool_root = pool;
    for (size_t i = 0; i < length_k; i++) {
        pool[i].right = &pool[i + 1];
        pool[i + 1].parent = &pool[i];
    }
    index = 0;
}

template <typename T>
void MaxMinSlideWindow<T>::add_new(const T& x) {
    MmswRBNode<T>*z, *y;
    if (index >= length_k) {
        y = DeleteNode(nodes[index % length_k]);
        pool_push(y);
    }
    z = pool_pop();
    z->key = x;
    z->index = static_cast<int>(index);
    z->size = 1;
    nodes[index % length_k] = z;
    InsertNode(z);
    index++;
}

#endif
