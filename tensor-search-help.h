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

float inline similarity_fun(const std::vector<float> &v1, const std::vector<float> &v2, const std::string &distance_type, const size_t v1_start = 0, size_t v1_count = 0)
{
    if (v1_count == 0)
    {
        v1_count = v1.size();
    }

    // if (v2_count == 0)
    // {
    //     v2_count = v2.size();
    // }

    // if (v1_count != v2_count)
    // {
    //     AU_EXIT("V1 and V2 have different size");
    // }

    // if (v1_start + v1_count > v1.size() || v2_start + v2_count > v2.size())
    // {
    //     std::cout << "v1_start = " << v1_start << "\n";
    //     std::cout << "v1_count = " << v1_count << "\n";
    //     std::cout << "v1.size() = " << v1.size() << "\n";

    //     std::cout << "v2_start = " << v2_start << "\n";
    //     std::cout << "v2_count = " << v2_count << "\n";
    //     std::cout << "v2.size() = " << v2.size() << "\n";

    //     AU_EXIT("Invalid start or count parameters");
    // }

    float result = 0.0;
    if (distance_type == "dotproduct")
    {
        for (size_t i = 0; i < v1_count; ++i)
        {
            result += v1[v1_start + i] * v2[i];
        }
    }
    else if (distance_type == "cosine")
    {
        float dot_product = 0.0;
        float norm_v1 = 0.0;
        float norm_v2 = 0.0;

        for (size_t i = 0; i < v1_count; ++i)
        {
            dot_product += v1[v1_start + i] * v2[i];
            norm_v1 += v1[v1_start + i] * v1[v1_start + i];
            norm_v2 += v2[i] * v2[i];
        }

        result = dot_product / (std::sqrt(norm_v1) * std::sqrt(norm_v2));
        result = 1.0 - result; // Cosine distance (1 - Cosine similarity)
    }
    else if (distance_type == "jaccard")
    {
        float intersection = 0.0;
        float union_size = 0.0;

        for (size_t i = 0; i < v1_count; ++i)
        {
            if (v1[v1_start + i] != 0 || v2[i] != 0)
            {
                union_size++;
                if (v1[v1_start + i] != 0 && v2[i] != 0)
                {
                    intersection++;
                }
            }
        }

        if (union_size == 0)
        {
            return 0.0; // Handle division by zero
        }

        result = 1.0 - (intersection / union_size); // Jaccard distance
    }
    else if (distance_type == "angular")
    {
        float dot_product = 0.0;
        float norm_v1 = 0.0;
        float norm_v2 = 0.0;

        for (size_t i = 0; i < v1_count; ++i)
        {
            dot_product += v1[v1_start + i] * v2[i];
            norm_v1 += v1[v1_start + i] * v1[v1_start + i];
            norm_v2 += v2[i] * v2[i];
        }

        if (norm_v1 == 0 || norm_v2 == 0)
        {
            return 0.0; // Handle division by zero
        }

        result = 1.0 - (dot_product / (std::sqrt(norm_v1) * std::sqrt(norm_v2))); // Angular distance
    }
    else if (distance_type == "euclidean")
    {
        for (size_t i = 0; i < v1_count; ++i)
        {
            float diff = v1[v1_start + i] - v2[i];
            result += diff * diff;
        }
        result = std::sqrt(result); // Euclidean distance
    }
    else
    {
        // Invalid distance type
        return 0.0;
    }

    return result;
}

void normalize_vector(std::vector<float> &data)
{
    // Calculate the magnitude (length) of the vector
    float magnitude = 0.0;
    for (float value : data)
    {
        magnitude += value * value;
    }
    magnitude = std::sqrt(magnitude);

    // Check for division by zero (if magnitude is very close to zero)
    if (magnitude < 1e-8)
    {
        // Handle this case as needed, e.g., by returning 0 or some default value
        // return 0.0;
        AU_EXIT("Can not normalize_vector");
    }

    // Normalize the vector by dividing each element by its magnitude
    float normalize_sum = 0;
    for (float &value : data)
    {
        value = value / magnitude;
    }
}

float average_vector(std::vector<float> &data)
{
    // Check if the vector is empty to avoid division by zero.
    if (data.empty())
    {
        return 0.0; // You can choose an appropriate value for an empty vector.
    }

    float sum = 0.0;
    for (const float &value : data)
    {
        sum += value;
    }

    return sum / static_cast<float>(data.size());
}

std::string inline extractDirFileName(const std::string &inputString)
{
    size_t colonPos = inputString.find(':');

    // Check if the colon is found
    if (colonPos != std::string::npos)
    {
        // Extract the substring before the colon
        return inputString.substr(0, colonPos);
    }
    else
    {
        // Return the entire input string if the colon is not found
        return inputString;
    }
}

bool inline isFile(const std::string &path)
{
    // Check if the path exists
    if (std::filesystem::exists(path))
    {
        // Check if the path is a regular file
        if (std::filesystem::is_regular_file(path))
        {
            return true; // If it's a file, return true
        }
        else
        {
            return false; // If it's a directory, return false
        }
    }
    else
    {
        std::cerr << "Path [" << path << "] does not exist." << std::endl;
        exit(1); // Exit the program with an error code
    }
}

template <typename T>
bool areVectorsEqual(const std::vector<T> &vector1, const std::vector<T> &vector2)
{
    if (vector1.size() != vector2.size())
    {
        return false; // Vectors have different sizes, so they can't be equal.
    }

    for (size_t i = 0; i < vector1.size(); i++)
    {
        if (vector1[i] != vector2[i])
        {
            return false; // Values at index i are not equal, so the vectors are not equal.
        }
    }

    return true; // All elements are equal, so the vectors are equal.
}

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
