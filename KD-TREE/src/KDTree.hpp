#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <size_t N, typename ElemType>
class KDTree {

 public:

  typedef std::pair<Point<N>, ElemType> value_type;

  KDTree();

  ~KDTree();

  KDTree(const KDTree &rhs);
  KDTree &operator=(const KDTree &rhs);

  size_t dimension() const;

  size_t size() const;
  bool empty() const;

  bool contains(const Point<N> &pt) const;

  void insert(const Point<N> &pt, const ElemType &value=ElemType());

  ElemType &operator[](const Point<N> &pt);

  ElemType &at(const Point<N> &pt);
  const ElemType &at(const Point<N> &pt) const;

private:

  struct KDNode{
    Point<N> pt_;
    KDNode* nodes[2];
    int level_;
    ElemType value_;

    KDNode(const Point<N>&pt, int level, const ElemType& value):pt_(pt),level_(level), value_(value) {nodes[0]=nodes[1]=nullptr;}
  };

  KDNode *root;
  size_t size_;

  KDNode *find(KDNode *curr_node, const Point<N>& pt) const;

  KDNode *copy(KDNode* root);

  void delete_recursive(KDNode* curr_node);

  
};

template <std::size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDNode* KDTree<N, ElemType>::find(typename KDTree<N, ElemType>::KDNode* curr_node, const Point<N>& pt) const {
    if (curr_node == NULL || curr_node->pt_ == pt) return curr_node;

    const Point<N>& currPoint = curr_node->pt_;
    int currLevel = curr_node->level_;
    if (pt[currLevel % N] < currPoint[currLevel % N]) {
        return curr_node->nodes[0] == NULL ? curr_node : find(curr_node->nodes[0], pt);
    } else {
        return curr_node->nodes[1] == NULL ? curr_node : find(curr_node->nodes[1], pt);
    }
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::KDNode* KDTree<N, ElemType>::copy(typename KDTree<N, ElemType>::KDNode* root_) {
    if (root_ == nullptr) return NULL;
    KDNode *new_root = new KDNode(*root_);
    new_root->nodes[0] = copy(new_root->nodes[0]);
    new_root->nodes[1] = copy(new_root->nodes[1]);
    return new_root;
}

template <std::size_t N, typename ElemType>
void KDTree<N, ElemType>::delete_recursive(typename KDTree<N, ElemType>::KDNode* curr_node) {
    if (curr_node == nullptr) return;
    delete_recursive(curr_node->nodes[0]);
    delete_recursive(curr_node->nodes[1]);
    --size_;
    delete curr_node;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() {
  this->root = nullptr; size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
  delete_recursive(this->root);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree &rhs) {
  this->root = copy(rhs.root);
  this->size_ = rhs.size_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType> &KDTree<N, ElemType>::operator=(const KDTree &rhs) {
  if (this != &rhs) {
        delete_recursive(this->root);
        this->root = copy(rhs.root);
        this->size_ = rhs.size_;
    }
  return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
  return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
  return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
  return size_==0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N> &pt) const {
  KDNode *node = find(this->root, pt);
  return node != nullptr && node->pt_ == pt;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    KDNode *targetNode = find(this->root, pt);
    if (targetNode == nullptr) { 
        this->root = new KDNode(pt, 0, value);
        this->size_ = 1;
    } else {
        if (targetNode->pt_ == pt) {
            targetNode->value_ = value;
        } else {
            int currLevel = targetNode->level_;
            KDNode* newNode = new KDNode(pt, currLevel + 1, value);
            if (pt[currLevel % N] < targetNode->pt_[currLevel % N]) {
                targetNode->nodes[0] = newNode;
            } else {
                targetNode->nodes[1] = newNode;
            }
            ++size_;
        }
    }
}

template <std::size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt) {
    KDNode *node = find(this->root, pt);
    if (node != nullptr && node->pt_ == pt) {
        return node->value_;
    } else {
        insert(pt);
        if (node == NULL) return this->root->value_;
        else return (node->nodes[0] != NULL && node->nodes[0]->pt_ == pt) ? node->nodes[0]->value_: node->nodes[1]->value_;
    }
}

template <std::size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    const KDTree<N, ElemType>& constThis = *this;
    return const_cast<ElemType&>(constThis.at(pt));
}

template <std::size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    KDNode *node = find(this->root, pt);
    if (node == nullptr || node->pt_ != pt) {
        throw std::out_of_range("Point not found in the KD-Tree");
    } else {
        return node->value_;
    }
}

#endif 
