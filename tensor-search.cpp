/**
 *
 * Author: Bin Dong dbin@lbl.gov
 * Web: https://crd.lbl.gov/bin-dong
 * Scientific Data Management Research Group
 * Lawrence Berkeley National Laboratory
 *
 */

#include <iostream>
#include <stdarg.h>
#include <vector>
#include <stdlib.h>
#include "ft.h"
#include "tensor-search-help.h"

using namespace std;
using namespace FT;

#define VEC_SIZE 8

std::string pattern_dir = "/Users/dbin/work/TensorSearch/db-pattern";
std::vector<std::vector<float>> pattern_data;
// Put all into config file
std::string A_endpoint_id;
std::string db_data_file_type = "EP_HDF5";
std::string db_data_dir = "/Users/dbin/work/TensorSearch/db-data";
std::string db_data_h5_dataset = "/testg/testd";
int input_data_type = 2; //  0 : short, 1 : double, 2 : float
std::string output_file_dir = "/Users/dbin/work/TensorSearch/tensor-search-result";
bool is_dir_output_match_replace_rgx = false;
std::string dir_output_match_rgx = "^(.*)\\.h5$";
std::string dir_output_replace_rgx = "$1-tensor-search.h5";
bool is_column_major = true;
std::vector<std::string> aug_merge_index;
int n_patterns;
std::vector<int> pattern_size;
int data_rank = 2;
std::vector<std::string> null_str;

// Put all into config file end

void all_gather_vector(const std::vector<float> &v_to_send, std::vector<float> &v_to_receive);

float inline similarity_fun(std::vector<float> &v1, std::vector<float> &v2)
{
    if (v1.size() != v2.size())
    {
        return 0.0;
    }

    float result = 0.0;
    for (size_t i = 0; i < v1.size(); ++i)
    {
        result += v1[i] * v2[i];
    }

    return result;
}

float normalize_similarity(std::vector<float> &data)
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
        return 0.0;
    }

    // Normalize the vector by dividing each element by its magnitude
    for (float &value : data)
    {
        value /= magnitude;
    }

    return magnitude; // Return the original magnitude for reference, if needed
}

template <class T1, class T2>
std::vector<float> calculate_similarity(std::vector<T1> &data, std::vector<int> &data_size, std::vector<T2> &pattern_data, std::vector<int> &pattern_data_size)
{
    // now we assume data_size == pattern_data_size
    int rank = data_size.size();
    assert(pattern_data_size.size() == rank);

    std::vector<float> similarity_result;
    similarity_result.resize(1, 0);
    if (rank == 1)
    {
        similarity_result[0] = similarity_fun(data, pattern_data);
    }
    else if (rank == 2)
    {
        size_t rows = data_size[0], cols = data_size[1];
        std::vector<std::vector<float>> data_2d = ConvertVector1DTo2D(data, rows, cols, false);
        std::vector<std::vector<float>> pattern_data_2d = ConvertVector1DTo2D(pattern_data, rows, cols, false);
        std::vector<float> similarity_result_temp;
        similarity_result_temp.resize(cols, 0);
        for (size_t i = 0; i < cols; i++)
        {
            similarity_result_temp[i] = similarity_fun(data_2d[i], pattern_data_2d[i]);
        }
        similarity_result[0] = normalize_similarity(similarity_result_temp);
    }
    else if (rank == 3)
    {
        std::vector<std::vector<std::vector<float>>> data_3d = ConvertVector1DTo3D(data, data_size[0], data_size[1], data_size[2], false);
        std::vector<std::vector<std::vector<float>>> pattern_data_3d = ConvertVector1DTo3D(pattern_data, pattern_data_size[0], pattern_data_size[1], pattern_data_size[2], false);
        std::vector<float> similarity_result_temp, similarity_result_temp2;
        size_t rows = data_size[0], cols = data_size[1], deps = data_size[2];
        similarity_result_temp.resize(cols, 0);
        similarity_result_temp2.resize(deps, 0);

        for (size_t j = 0; j < deps; j++)
        {
            for (size_t i = 0; i < cols; i++)
            {
                similarity_result_temp[i] = similarity_fun(data_3d[j][i], pattern_data_3d[j][i]);
            }
            similarity_result_temp2[j] = normalize_similarity(similarity_result_temp);
        }
        similarity_result[0] = normalize_similarity(similarity_result_temp2);
    }
    else
    {
        AU_EXIT("We haven't implement this 4D and beyond yet");
    }

    return similarity_result;
}

template <class TT>
inline Stencil<std::vector<float>> tensor_search_udf(const Stencil<TT> &iStencil)
{
    std::vector<int> max_offset_upper; // Size of input data for the whole chunk, zero based
    iStencil.GetOffsetUpper(max_offset_upper);

    if (ft_rank == 0 || ft_rank == (ft_size - 1))
        PrintVector("max_offset_upper = ", max_offset_upper);
    std::vector<int> start_offset;
    start_offset.resize(max_offset_upper.size(), 0);

    std::vector<TT> db_data_per_udf;
    iStencil.ReadNeighbors(start_offset, max_offset_upper, db_data_per_udf);

    // std::vector<std::vector<float>> ts2d;
    // ts2d = DasLib::Vector1D2D(chs_per_file_udf, db_data_per_udf);

    Stencil<std::vector<float>> oStencil;
    std::vector<std::vector<float>> similarity_vectors;
    std::vector<size_t> vector_shape(2);

    similarity_vectors.resize(n_patterns);
    for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
    {
        similarity_vectors[pattern_index] = calculate_similarity(db_data_per_udf, max_offset_upper, pattern_data[pattern_index], pattern_size);
    }
    vector_shape[0] = n_patterns;
    vector_shape[1] = similarity_vectors[0].size();

    std::vector<float> similarity_vector_1d = ConvertVector2DTo1D(similarity_vectors);

    oStencil.SetShape(vector_shape);
    oStencil = similarity_vector_1d;
    return oStencil;
}

void load_process_pattern_files_on_each_process()
{
    AU::Array<float> *pattern_arrays;
    std::vector<int> chunk_size, overlap_size = {0, 0};
    std::vector<std::string> file_size_str_h5;

    pattern_arrays = new AU::Array<float>("EP_DIR:EP_HDF5:" + pattern_dir + ":/testg/testd");
    pattern_arrays->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);
    pattern_arrays->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str_h5);
    String2Vector(file_size_str_h5[0], chunk_size);
    assert(chunk_size.size() == data_rank);
    overlap_size.resize(data_rank, 0);
    pattern_arrays->SetChunkSize(chunk_size);
    pattern_arrays->SetOverlapSize(overlap_size);
    pattern_size = chunk_size;

    if (!ft_rank)
    {
        PrintVector("db-pattern: chunk_size = ", chunk_size);
        PrintVector("db-pattern: overlap_size = ", overlap_size);
    }

    int n_patterns_on_my_rank = 0;
    if (ft_size > 1)
    {
        pattern_arrays->SetChunkSchedulingMethod(CHUNK_SCHEDULING_CR);
        unsigned long long my_chunk_start, my_chunk_end;
        pattern_arrays->GetMyChunkStartEnd(my_chunk_start, my_chunk_end);
        if (!ft_rank)
            std::cout << "rank [" << ft_rank << "]: my_chunk_start = " << my_chunk_start << ", my_chunk_end = " << my_chunk_end << "\n";
        n_patterns_on_my_rank = my_chunk_end - my_chunk_start;
    }
    pattern_arrays->ControlEndpoint(DIR_N_FILES, file_size_str_h5);
    n_patterns = std::stoi(file_size_str_h5[0]);

    int n_patterns_to_go;
    if (ft_size > 1)
    {
        n_patterns_to_go = n_patterns_on_my_rank;
    }
    else
    {
        n_patterns_to_go = n_patterns;
    }

    if (!ft_rank)
        std::cout << "Rank " << ft_rank << "# of total pattern files = " << n_patterns << ", # of pattern files per process = " << n_patterns_to_go << std::endl;

    pattern_data.resize(n_patterns_to_go);

    std::vector<float> pattern_data_per_file;
    for (int pattern_file_index = 0; pattern_file_index < n_patterns_to_go; pattern_file_index++)
    {
        pattern_arrays->ReadNextChunk(pattern_data_per_file);

        pattern_data[pattern_file_index] = pattern_data_per_file;
    }

    if (ft_size > 1)
    {
        int points_per_pattern, max_points_per_pattern;
        if ((!pattern_data.empty()) && (!pattern_data[0].empty()))
        {
            points_per_pattern = pattern_data[0].size();
        }

        MPI_Allreduce(&points_per_pattern, &max_points_per_pattern, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        // MPI_Allreduce(&timepoints, &max_timepoints, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

        std::vector<float> pattern_data_1D = ConvertVector2DTo1D(pattern_data);
        std::vector<float> pattern_data_1D_merged;
        all_gather_vector(pattern_data_1D, pattern_data_1D_merged);
        pattern_data = ConvertVector1DTo2D(pattern_data_1D_merged, n_patterns, max_points_per_pattern, true);
    }
}
int main(int argc, char *argv[])
{
    // Init the MPICH, etc.
    AU_Init(argc, argv);

    A_endpoint_id = "EP_DIR:" + db_data_file_type + ":" + db_data_dir + ":" + db_data_h5_dataset;

    // Input data
    // AU::Array<short> *

    FT::ArrayBase *A;
    AuEndpointDataType t = AuEndpointDataType::AU_FLOAT;
    if (input_data_type == 0)
    {
        t = AuEndpointDataType::AU_SHORT;
        A = new FT::Array<short>(A_endpoint_id);
    }
    else if (input_data_type == 1)
    {
        t = AuEndpointDataType::AU_DOUBLE;
        A = new FT::Array<double>(A_endpoint_id);
    }
    else if (input_data_type == 2)
    {
        t = AuEndpointDataType::AU_FLOAT;
        A = new FT::Array<float>(A_endpoint_id);
    }
    else
    {
        std::cout << "Not supported input_data_type \n";
        exit(-1);
    }

    A->SkipFileTail();
    A->ExecuteUDFOnce();
    A->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);

    std::vector<int> chunk_size;
    std::vector<int> overlap_size;

    // if (is_column_major)
    // {
    //     aug_merge_index.push_back("0");
    // }
    // else
    // {
    //     aug_merge_index.push_back("1");
    // }
    // A->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);

    std::vector<std::string> file_size_str;
    A->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
    String2Vector(file_size_str[0], chunk_size);

    data_rank = chunk_size.size();
    overlap_size.resize(data_rank, 0);
    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);

    if (!ft_rank)
    {
        PrintVector("chunk_size = ", chunk_size);
        PrintVector("overlap_size = ", overlap_size);
    }

    AU::Array<float> *B;
    B = new AU::Array<float>("EP_DIR:" + db_data_file_type + ":" + output_file_dir + ":" + db_data_h5_dataset);
    // Use the below rgx pattern to name the file
    std::vector<std::string> aug_output_replace_arg;
    aug_output_replace_arg.push_back(dir_output_match_rgx);
    aug_output_replace_arg.push_back(dir_output_replace_rgx);
    // B->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);

    if (is_dir_output_match_replace_rgx)
        B->ControlEndpoint(DIR_OUPUT_REPLACE_RGX, aug_output_replace_arg);

    A->EnableApplyStride(chunk_size);
    A->SetVectorDirection(AU_FLAT_OUTPUT_ROW);

    load_process_pattern_files_on_each_process();

    // TRANSFORM_NO_MP(A, tensor_search_udf, B, t, std::vector<float>);

    static_cast<FT::Array<float> *>(A)->TransformNoMP<std::vector<float>>(tensor_search_udf, B);

    A->ReportCost();

    // Clear
    delete A;
    delete B;

    AU_Finalize();

    return 0;
}

void all_gather_vector(const std::vector<float> &v_to_send, std::vector<float> &v_to_receive)
{
    // https://rookiehpc.github.io/mpi/docs/mpi_allgatherv/index.html
    int local_sum = v_to_send.size();

    std::vector<int> local_sum_vector;
    local_sum_vector.resize(ft_size);
    MPI_Allgather(&local_sum, 1, MPI_INT, local_sum_vector.data(), 1, MPI_INT, MPI_COMM_WORLD);

    int global_sum = 0;
    std::vector<int> displacements;
    displacements.resize(ft_size);
    displacements[0] = 0;
    for (int i = 0; i < local_sum_vector.size(); i++)
    {
        global_sum = global_sum + local_sum_vector[i];
        if (i > 0)
            displacements[i] = displacements[i - 1] + local_sum_vector[i - 1];
    }
    v_to_receive.resize(global_sum);
    MPI_Allgatherv(v_to_send.data(), local_sum, MPI_FLOAT, v_to_receive.data(), local_sum_vector.data(), displacements.data(), MPI_FLOAT, MPI_COMM_WORLD);
}
