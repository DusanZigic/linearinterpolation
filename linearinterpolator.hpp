#ifndef LINEAR_INTERPOLATOR_HPP
#define LINEAR_INTERPOLATOR_HPP

#include <vector>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <cmath>

template<typename T>
class LinearInterpolator {
public:
	LinearInterpolator();
	LinearInterpolator(const std::vector<T> &xData, const std::vector<T> &fData);
	LinearInterpolator(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &fData);
	LinearInterpolator(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<std::vector<T>> &fData);
	LinearInterpolator(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &fData);
	LinearInterpolator(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &x4Data, const std::vector<T> &fData);
	
	void setData(const std::vector<T> &xData, const std::vector<T> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<std::vector<T>> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &fData);
	void setData(const std::vector<T> &x1Data, const std::vector<T> &x2Data, const std::vector<T> &x3Data, const std::vector<T> &x4Data, const std::vector<T> &fData);

private:
	std::vector<std::vector<T>> m_axes;
    std::vector<T> m_values;
    std::array<size_t, 3> m_strides;
    int m_dim;

	void precomputeStrides();
};

#endif //LINEAR_INTERPOLATOR_HPP