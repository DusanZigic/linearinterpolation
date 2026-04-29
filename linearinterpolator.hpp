#ifndef LINEAR_INTERPOLATOR_HPP
#define LINEAR_INTERPOLATOR_HPP

#include <vector>
#include <array>
#include <algorithm>
// #include <stdexcept>

template<typename T>
class LinearInterpolator {
public:
	LinearInterpolator(std::string name = "unnamed") : m_name(name), m_dim(0) {}

	LinearInterpolator(const std::vector<T> &x, const std::vector<T> &f, std::string name = "unnamed") : m_name(name) {
		initializeInternal({x}, f);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &f, std::string name = "unnamed") : m_name(name) { 
		initializeInternal({x1, x2}, f);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<std::vector<T>> &f_grid, std::string name = "unnamed") : m_name(name) {
		std::vector<T> f_flat;
        f_flat.reserve(x1.size() * x2.size());
        for (const auto& row : f_grid) {
            f_flat.insert(f_flat.end(), row.begin(), row.end());
        }
        initializeInternal({x1, x2}, f_flat);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &f, std::string name = "unnamed") : m_name(name) { 
		initializeInternal({x1, x2, x3}, f);
	}

	LinearInterpolator(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &x4, const std::vector<T> &f, std::string name = "unnamed") : m_name(name) { 
		initializeInternal({x1, x2, x3, x4}, f);
	}

	void setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &f) {
		initializeInternal({x1, x2}, f);
	}

    void setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &f) {
		initializeInternal({x1, x2, x3}, f);
	}

    void setData(const std::vector<T> &x1, const std::vector<T> &x2, const std::vector<T> &x3, const std::vector<T> &x4, const std::vector<T> &f) {
		initializeInternal({x1, x2, x3, x4}, f);
	}

	inline T interpolate(T x) const {
        size_t i = locateGridIndex(0, x);
        size_t ss = m_searchStrides[0];
        
        T x0 = m_axes[0][i * ss], x1 = m_axes[0][(i + 1) * ss];
        T frac = (x - x0) / (x1 - x0);
        return m_values[i] * (1 - frac) + m_values[i + 1] * frac;
    }

    inline T interpolate(T x1, T x2) const {
        size_t i1 = locateGridIndex(0, x1);
        size_t i2 = locateGridIndex(1, x2);

        size_t ss1 = m_searchStrides[0], ss2 = m_searchStrides[1];
        size_t s0 = m_memStrides[0];

        T w1 = (x1 - m_axes[0][i1 * ss1]) / (m_axes[0][(i1 + 1) * ss1] - m_axes[0][i1 * ss1]);
        T w2 = (x2 - m_axes[1][i2 * ss2]) / (m_axes[1][(i2 + 1) * ss2] - m_axes[1][i2 * ss2]);

        size_t idx = i1 * s0 + i2;
        T v00 = m_values[idx];
        T v01 = m_values[idx + 1];
        T v10 = m_values[idx + s0];
        T v11 = m_values[idx + s0 + 1];

        return (1-w1)*(1-w2)*v00 + (1-w1)*w2*v01 + w1*(1-w2)*v10 + w1*w2*v11;
    }

    inline T interpolate(T x1, T x2, T x3) const {
        size_t i1 = locateGridIndex(0, x1), i2 = locateGridIndex(1, x2), i3 = locateGridIndex(2, x3);
        size_t ss1 = m_searchStrides[0], ss2 = m_searchStrides[1], ss3 = m_searchStrides[2];
        size_t s0 = m_memStrides[0], s1 = m_memStrides[1];

        T w1 = (x1 - m_axes[0][i1*ss1]) / (m_axes[0][(i1+1)*ss1] - m_axes[0][i1*ss1]);
        T w2 = (x2 - m_axes[1][i2*ss2]) / (m_axes[1][(i2+1)*ss2] - m_axes[1][i2*ss2]);
        T w3 = (x3 - m_axes[2][i3*ss3]) / (m_axes[2][(i3+1)*ss3] - m_axes[2][i3*ss3]);

        auto g = [&](size_t o1, size_t o2, size_t o3) { 
            return m_values[(i1+o1)*s0 + (i2+o2)*s1 + (i3+o3)]; 
        };

        T c00 = g(0,0,0)*(1-w1) + g(1,0,0)*w1;
        T c01 = g(0,0,1)*(1-w1) + g(1,0,1)*w1;
        T c10 = g(0,1,0)*(1-w1) + g(1,1,0)*w1;
        T c11 = g(0,1,1)*(1-w1) + g(1,1,1)*w1;

        T c0 = c00*(1-w2) + c10*w2;
        T c1 = c01*(1-w2) + c11*w2;

        return c0*(1-w3) + c1*w3;
    }

    inline T interpolate(T x1, T x2, T x3, T x4) const {
        size_t i[4] = {locateGridIndex(0, x1), locateGridIndex(1, x2), locateGridIndex(2, x3), locateGridIndex(3, x4)};
        T w[4];
        for(int d=0; d<4; ++d) {
            size_t base = i[d] * m_searchStrides[d];
            size_t next = (i[d] + 1) * m_searchStrides[d];
            w[d] = ( (d==0?x1:(d==1?x2:(d==2?x3:x4))) - m_axes[d][base] ) / (m_axes[d][next] - m_axes[d][base]);
        }

        T res = 0;
        for (int b = 0; b < 16; ++b) {
            T wt = 1.0;
            size_t idx = 0;
            for (int d = 0; d < 4; ++d) {
                bool bit = (b >> (3 - d)) & 1;
                wt *= bit ? w[d] : (1.0 - w[d]);
                idx += (i[d] + (bit ? 1 : 0)) * (d == 3 ? 1 : m_memStrides[d]);
            }
            res += wt * m_values[idx];
        }
        return res;
    }

private:
	std::string m_name;
	std::vector<std::vector<T>> m_axes;
    std::vector<T> m_values;
    std::vector<size_t> m_gridSizes;
    std::vector<size_t> m_memStrides;
    std::vector<size_t> m_searchStrides;
    int m_dim;
	
	void initializeInternal(const std::vector<std::vector<T>> &axes, const std::vector<T> &f) {
        m_axes = axes;
        m_values = f;
        m_dim = static_cast<int>(axes.size());
        precomputeStrides();
    }

	void precomputeStrides() {
		m_gridSizes.resize(m_dim);
        m_searchStrides.resize(m_dim);
        m_memStrides.assign(m_dim, 1);

        for (int d = 0; d < m_dim; ++d) {
			const auto& vec = m_axes[d];

            size_t step = 1;
            while (step < m_axes[d].size() && m_axes[d][step] == m_axes[d][0]) step++;
            m_searchStrides[d] = step;

			if (vec.size() == m_values.size()) {
				size_t period = vec.size();
				for (size_t i = step; i < vec.size(); ++i) {
					if (vec[i] == vec[0] && vec[i-1] != vec[0]) {
						period = i;
						break;
					}
				}
				m_gridSizes[d] = period / step;
			} else {
				m_searchStrides[d] = 1;
                m_gridSizes[d] = vec.size();
			}
        }

        for (int d = m_dim - 2; d >= 0; --d) {
            m_memStrides[d] = m_gridSizes[d + 1] * m_memStrides[d + 1];
        }
	}

	inline size_t locateGridIndex(int axis, T val) const {
        const auto& vec = m_axes[axis];
        size_t ss = m_searchStrides[axis];
        size_t n = m_gridSizes[axis];

        if (val <= vec[0]) return 0;
        if (val >= vec[(n - 1) * ss]) return n - 2;

        size_t low = 0, high = n - 2, ans = 0;
        while (low <= high) {
            size_t mid = low + (high - low) / 2;
            if (vec[mid * ss] <= val) {
                ans = mid;
                low = mid + 1;
            } else {
                high = mid - 1;
            }
        }
        return ans;
    }
};

#endif //LINEAR_INTERPOLATOR_HPP