struct SimilarityStruct
{
    AU_UDT_INIT(SimilarityStruct)
    float similarity;
    unsigned long long index;
};

void PrintSSV(const std::string &note, const std::vector<SimilarityStruct> &ssv)
{
    std::cout << "Index : Similarity of " << note << " at " << ft_rank << "\n";
    for (int i = 0; i < ssv.size(); i++)
    {
        std::cout << ssv[i].index << ": " << ssv[i].similarity << "\n";
    }
}

bool compareSimilarity(const SimilarityStruct &a, const SimilarityStruct &b)
{
    return a.similarity > b.similarity;
}

std::vector<SimilarityStruct> FindTopKAfterCollect(const std::vector<SimilarityStruct> &v, int k, int n_patterns)
{
    std::vector<SimilarityStruct> result_final;

    result_final.reserve(n_patterns * k);
    // Ensure k is not greater than the size of the input vector
    k = std::min(k, static_cast<int>(v.size()));

    for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
    {
        std::vector<SimilarityStruct> result;
        // Iterate over segments
        for (int i = pattern_index * k; i <= v.size() - k; i += (n_patterns * k))
        {
            result.insert(result.end(), v.begin() + i, v.begin() + i + k);
        }

        // Sort the result vector in descending order to get the overall top k elements
        std::sort(result.begin(), result.end(), compareSimilarity);
        // PrintSSV("template " + std::to_string(pattern_index), result);
        if (result.size() > k)
        {
            result.resize(k);
        }

        // Resize the result vector to contain only the top k elements
        // result.resize(k);
        result_final.insert(result_final.end(), result.begin(), result.end());
        // PrintSSV("template " + std::to_string(pattern_index), result_final);
    }
    return result_final;
}

std::vector<SimilarityStruct> CollectSimilarityStructToRoot(const std::vector<SimilarityStruct> &v)
{
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int num_elements_per_rank = v.size();
    int root = 0;

    // Determine the total number of elements
    int total_elements = size * num_elements_per_rank;

    // Allocate memory for the receive buffer on the root process
    std::vector<SimilarityStruct> gathered_data(total_elements);

    // Perform the gather operation
    MPI_Gather(v.data(), num_elements_per_rank * sizeof(SimilarityStruct), MPI_BYTE,
               gathered_data.data(), num_elements_per_rank * sizeof(SimilarityStruct), MPI_BYTE,
               root, MPI_COMM_WORLD);

    if (rank == root)
    {
        return gathered_data;
    }
    else
    {
        return {}; // Non-root processes return an empty vector
    }
}

// Function to extract the top K largest elements and their indices
std::vector<SimilarityStruct> top_k_largest(const std::vector<float> &v, size_t k, unsigned long long my_index_start)
{
    // Ensure the input vector is not empty
    if (v.empty() || k == 0)
    {
        std::cerr << "Invalid input parameters." << std::endl;
        return {};
    }

    // Create a vector of SimilarityStruct to store values along with their indices
    std::vector<SimilarityStruct> similarityStructs;
    for (unsigned long long i = 0; i < v.size(); ++i)
    {
        similarityStructs.push_back({v[i], my_index_start + i});
    }

    // // Sort the vector of SimilarityStruct in descending order based on similarity
    // std::sort(similarityStructs.begin(), similarityStructs.end(),
    //           [](const SimilarityStruct &a, const SimilarityStruct &b)
    //           {
    //               return a.similarity > b.similarity;
    //           });

    // Use partial_sort to get the top K elements
    std::partial_sort(similarityStructs.begin(), similarityStructs.begin() + k, similarityStructs.end(),
                      [](const SimilarityStruct &a, const SimilarityStruct &b)
                      {
                          return a.similarity > b.similarity;
                      });
    // Resize the result vector to hold the top K elements
    std::vector<SimilarityStruct> topKElements(k);

    // Copy the top K elements
    std::copy(similarityStructs.begin(), similarityStructs.begin() + k, topKElements.begin());

    return topKElements;
}

// Function to insert a new element into the top-k vector
void insert_topk(std::vector<SimilarityStruct> &v, const float similarity, const unsigned long long index)
{
    // Find the smallest similarity in the current top-k vector
    float smallestSimilarity = v.back().similarity;

    // If the new similarity is greater than the smallest, insert the new element
    if (similarity > smallestSimilarity)
    {
        // Find the position to insert the new element while maintaining sorted order
        auto it = std::lower_bound(v.begin(), v.end(), similarity,
                                   [](const SimilarityStruct &s1, float value)
                                   {
                                       return s1.similarity > value;
                                   });

        // Insert the new element at the found position
        v.insert(it, {similarity, index});

        // Remove the last element to keep the vector size as top_k
        v.pop_back();
    }
}

template <class T>
std::vector<std::vector<T>> Reshape2DVector(const std::vector<std::vector<T>> &nums, int m, int n)
{
    // Get the number of rows and columns in the original matrix
    int originalRows = nums.size();
    int originalCols = (originalRows > 0) ? nums[0].size() : 0;

    // Check if reshape is possible
    if (originalRows * originalCols != m * n)
    {
        std::cout << "Reshape is not possible.\n";
        return nums; // Return the original matrix if reshape is not possible
    }

    // Create the reshaped matrix with dimensions m x n
    std::vector<std::vector<T>> reshaped(m, std::vector<T>(n, T()));

    // Populate the reshaped matrix with elements from the original matrix using row-major order
    for (int i = 0; i < m * n; ++i)
    {
        reshaped[i / n][i % n] = nums[i / originalCols][i % originalCols];
    }

    return reshaped;
}

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

/**
 * @brief Convert a N x M matrix src to M x N matrix dst
 *
 * @tparam T
 * @param src
 * @param dst
 * @param N : row
 * @param M : column
 */

template <class T>
inline void transpose(T *src, T *dst, const int N, const int M)
{
    int i, j;
    for (int n = 0; n < N * M; n++)
    {
        i = n / N;
        j = n % N;
        dst[n] = src[M * j + i];
    }
}

template <class T>
inline std::vector<T> ConvertVector2DTo1D(const std::vector<std::vector<T>> &data2d)
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

template <class T1>
inline std::vector<SimilarityStruct> ConvertVector2DTo1DTopK(const std::vector<std::vector<T1>> &data2d, const size_t top_k, unsigned long long my_index_start)
{
    if (top_k == 0)
    {
        std::vector<SimilarityStruct> result;
        size_t rows = data2d.size();
        size_t cols = data2d[0].size();

        result.reserve(rows * cols);
        for (std::size_t i = 0; i < rows; ++i)
        {
            for (std::size_t j = 0; j < cols; ++j)
            {
                result.push_back({static_cast<float>(data2d[i][j]), j + my_index_start});
            }

            // copy(data2d[i].begin(), data2d[i].end(), back_inserter(result));
        }
        assert(result.size() == (rows * cols));
        return result;
    }
    else
    {

        std::vector<SimilarityStruct> result;
        size_t rows = data2d.size();
        result.reserve(rows * top_k);
        std::vector<SimilarityStruct> result_topk_per_row;
        for (std::size_t i = 0; i < rows; ++i)
        {
            result_topk_per_row = top_k_largest(data2d[i], top_k, my_index_start);
            copy(result_topk_per_row.begin(), result_topk_per_row.end(), back_inserter(result));
        }
        assert(result.size() == (rows * top_k));
        return result;
    }
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
        // for (size_t i = 0; i < v1_count; ++i)
        // {
        //     result += v1[v1_start + i] * v2[i];
        // }

        size_t i;
        for (i = 0; i < v1_count && i + 3 < v1_count; i += 4)
        {
            result += v1[v1_start + i] * v2[i] +
                      v1[v1_start + i + 1] * v2[i + 1] +
                      v1[v1_start + i + 2] * v2[i + 2] +
                      v1[v1_start + i + 3] * v2[i + 3];
        }

        // Handle remaining elements
        for (; i < v1_count; ++i)
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
