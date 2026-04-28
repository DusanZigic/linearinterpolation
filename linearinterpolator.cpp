#include "linearinterpolator.hpp"

template <typename T>
LinearInterpolator<T>::LinearInterpolator() : m_dim(0) {}

template <typename T>
LinearInterpolator<T>::LinearInterpolator(const std::vector<T> &x, const std::vector<T> &f) {
    setData(x, f);
}

template <typename T>
LinearInterpolator<T>::LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &f) {
    setData(x1, x2, f);
}

template <typename T>
LinearInterpolator<T>::LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<std::vector<T>> &f) {
    setData(x1, x2, f);
}

template <typename T>
LinearInterpolator<T>::LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &f) {
    setData(x1, x2, x3, f);
}

template <typename T>
LinearInterpolator<T>::LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &x4, const std::vector<T> &f) {
    setData(x1, x2, x3, x4, f);
}

template <typename T>
void LinearInterpolator<T>::setData(const std::vector<T> &x, const std::vector<T> &f) {
    m_axes = {x};
    m_values = f;
    m_dim = 1;
    precomputeStrides();
}

template <typename T>
void LinearInterpolator<T>::setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &f) {
    m_axes = {x1, x2};
    m_values = f;
    m_dim = 2;
    precomputeStrides();
}

template <typename T>
void LinearInterpolator<T>::setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<std::vector<T>> &f) {
    m_axes = {x1, x2};
    m_values.clear();
    m_values.reserve(x1.size() * x2.size());
    for(const auto& row : f) {
        m_values.insert(m_values.end(), row.begin(), row.end());
    }
    m_dim = 2;
    precomputeStrides();
}

template <typename T>
void LinearInterpolator<T>::setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &f) {
    m_axes = {x1, x2, x3};
    m_values = f;
    m_dim = 3;
    precomputeStrides();
}

template <typename T>
void LinearInterpolator<T>::setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &x4, const std::vector<T> &f) {
    m_axes = {x1, x2, x3, x4};
    m_values = f;
    m_dim = 4;
    precomputeStrides();
}

template <typename T>
void LinearInterpolator<T>::precomputeStrides() {
    if (m_dim < 2) return;

    m_strides[m_dim-2] = m_axes[m_dim-1].size();
    for (int i = m_dim - 3; i >= 0; --i) 
        m_strides[i] = m_axes[i+1].size() * m_strides[i+1];
}