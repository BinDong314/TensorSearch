template <class T>
inline std::vector<std::vector<T>> ConvertVector1DTo2D(const std::vector<T> &data1d, const size_t rows, const size_t cols, const bool row_layout_flag)
{
    std::vector<std::vector<T>> result;
    if (row_layout_flag)
    {
        // size_t rows = data1d.size() / cols;
        result.resize(rows);
        for (std::size_t i = 0; i < rows; ++i)
        {
            result[i].resize(cols);
            for (std::size_t j = 0; j < cols; ++j)
            {
                result[i][j] = data1d[i * cols + j];
            }
        }
    }
    else
    {
        // size_t rows = data1d.size() / cols;
        result.resize(cols);
        for (std::size_t i = 0; i < cols; ++i)
        {
            result[i].resize(rows);
            for (std::size_t j = 0; j < rows; ++j)
            {
                result[i][j] = data1d[j * cols + i];
            }
        }
    }

    return result;
}

template <class T>
inline std::vector<T> ConvertVector2DTo1D(std::vector<std::vector<T>> &data2d)
{
    std::vector<T> result;
    size_t rows = data2d.size();
    size_t cols = data2d[0].size();

    result.reserve(rows * cols);
    for (std::size_t i = 0; i < rows; ++i)
    {
        copy(data2d[i].begin(), data2d[i].end(), back_inserter(result));
    }
    assert(result.size() == (rows * cols));
    return result;
}

#include <vector>
#include <cassert>

template <class T>
inline std::vector<std::vector<std::vector<T>>> ConvertVector1DTo3D(const std::vector<T> &data1d, const size_t dim1, const size_t dim2, const size_t dim3, const bool row_layout_flag)
{
    std::vector<std::vector<std::vector<T>>> result;
    if (row_layout_flag)
    {
        // Calculate the number of elements in each dimension.
        size_t rows = dim1;
        size_t cols = dim2;
        size_t depth = dim3;

        result.resize(rows);
        for (std::size_t i = 0; i < rows; ++i)
        {
            result[i].resize(cols);
            for (std::size_t j = 0; j < cols; ++j)
            {
                result[i][j].resize(depth);
                for (std::size_t k = 0; k < depth; ++k)
                {
                    result[i][j][k] = data1d[i * (cols * depth) + j * depth + k];
                }
            }
        }
    }
    else
    {
        // Calculate the number of elements in each dimension.
        size_t rows = dim1;
        size_t cols = dim2;
        size_t depth = dim3;

        result.resize(depth);
        for (std::size_t k = 0; k < depth; ++k)
        {
            result[k].resize(cols);
            for (std::size_t j = 0; j < cols; ++j)
            {
                result[k][j].resize(rows);
                for (std::size_t i = 0; i < rows; ++i)
                {
                    result[k][j][i] = data1d[i * (cols * depth) + j * depth + k];
                }
            }
        }
    }

    return result;
}

template <class T>
inline std::vector<T> ConvertVector3DTo1D(std::vector<std::vector<std::vector<T>>> &data3d)
{
    std::vector<T> result;
    size_t dim1 = data3d.size();
    size_t dim2 = data3d[0].size();
    size_t dim3 = data3d[0][0].size();

    result.reserve(dim1 * dim2 * dim3);
    for (std::size_t i = 0; i < dim1; ++i)
    {
        for (std::size_t j = 0; j < dim2; ++j)
        {
            for (std::size_t k = 0; k < dim3; ++k)
            {
                result.push_back(data3d[i][j][k]);
            }
        }
    }
    assert(result.size() == (dim1 * dim2 * dim3));
    return result;
}
