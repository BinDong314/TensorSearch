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
#include <filesystem>

#include "tensor-search-help.h"

using namespace std;
using namespace FT;

#define VEC_SIZE 8

std::string pattern_dir = "/Users/dbin/work/TensorSearch/db-pattern-2d::/testg/testd";

std::vector<std::vector<float>> pattern_data;
// Put all into config file
std::string A_endpoint_id;
std::string db_data_file_type = "EP_HDF5";
std::string db_data_dir = "/Users/dbin/work/TensorSearch/db-data-2d:/testg/testd";
// std::string db_data_h5_dataset = "/testg/testd";
int input_data_type = 2; //  0 : short, 1 : double, 2 : float
std::string output_file_dir = "/Users/dbin/work/TensorSearch/db-data-pattern-2d:/testg/testd";
bool is_dir_output_match_replace_rgx = false;
std::string dir_output_match_rgx = "^(.*)\\.h5$";
std::string dir_output_replace_rgx = "$1-tensor-search.h5";
bool is_column_major = true;
std::vector<std::string> aug_merge_index;
int n_patterns;
std::vector<int> each_pattern_size;
int data_rank = 2;
int user_similarity_rank = 2;
bool user_has_similarity_rank = false;
std::vector<std::string> null_str;
std::string distance_type = "dotproduct";
bool is_input_single_file = false;
bool is_shift = false;
int shift_size = 1;
// Put all into config file end

void all_gather_vector(const std::vector<float> &v_to_send, std::vector<float> &v_to_receive);
void printf_help(char *cmd);

template <class T>
void calculate_similarity(const std::vector<T> &data, const std::vector<int> &data_size, const std::vector<std::vector<T>> &pattern_data, const std::vector<int> &each_pattern_data_size, const int n_patterns, std::vector<std::vector<float>> &similarity_result)
{
    /// similarity_result.resize(n_patterns);

    similarity_result.resize(n_patterns, std::vector<float>(1, 0.0));

    // now we assume data_size == pattern_data_size
    int final_similarity_rank;

    if (user_has_similarity_rank && user_similarity_rank == 1 && data_size.size() == 2)
    {
        final_similarity_rank = user_similarity_rank;
    }
    else
    {
        final_similarity_rank = data_size.size();
    }
    // assert(each_pattern_data_size.size() == final_similarity_rank);

    // if (!areVectorsEqual(each_pattern_data_size, data_size))
    // {
    //     PrintVector("each_pattern_data_size=", each_pattern_data_size);
    //     PrintVector("data_size = ", data_size);
    //     AU_EXIT("calculate_similarity has different input size!");
    // }

    // std::vector<float> similarity_result_of_each_pattern;
    // similarity_result_of_each_pattern.resize(1, 0);

    if (final_similarity_rank == 1)
    {
        // std::vector<float> data_1d = data;
        // std::vector<float> pattern_data_1d = pattern_data[pattern_index];
        // PrintVector("data = ", data_1d);
        // PrintVector("pattern_data = ", pattern_data_1d);
        if (!is_shift)
        {
            for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
            {
                similarity_result[pattern_index][0] = similarity_fun(data, pattern_data[pattern_index], distance_type);
                PrintVector("similarity_result = ", similarity_result[pattern_index]);
            }
        }
        else
        {
            size_t final_shift_size;
            if (shift_size == -1)
            {
                final_shift_size = pattern_data[0].size();
            }
            else
            {
                final_shift_size = shift_size;
            }
            int data_n_rows = data_size[0];
            // int n_pairs_of_vector = data_n_rows;
            std::cout << "data.size() = " << data.size() << "\n";
            std::cout << "data_n_rows = " << data_n_rows << "\n";
            std::cout << "final_shift_size = " << final_shift_size << "\n";
            // std::cout << "n_pairs_of_vector = " << data_n_rows << "\n";
            // std::size_t data_start_offset = 0;
            for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
            {
                similarity_result[pattern_index].resize(data_n_rows);
            }
            for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
            {
                if ((pattern_index % 10 == 0))
                    std::cout << "pattern_index = " << pattern_index << "\n";
// data_start_offset = 0;
#if defined(_OPENMP)
#pragma omp parallel for
#endif
                for (int shift_index = 0; shift_index < data_n_rows; shift_index++)
                {
                    similarity_result[pattern_index][shift_index] = similarity_fun(data, pattern_data[pattern_index], distance_type, shift_index * final_shift_size, final_shift_size);
                    // similarity_fun(data, pattern_data[pattern_index], distance_type, shift_index * final_shift_size, final_shift_size);
                    //  data_start_offset = data_start_offset + final_shift_size;
                }
            }
        }
    }
    else if (final_similarity_rank == 2)
    {
        size_t rows = data_size[0], cols = data_size[1];
        std::vector<std::vector<float>> data_2d = ConvertVector1DTo2D(data, rows, cols, false);

        if (!is_shift)
        {
            for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
            {
                std::vector<std::vector<float>> pattern_data_2d = ConvertVector1DTo2D(pattern_data[pattern_index], rows, cols, false);
                PrintVector("data = ", data);
                PrintVV("data_2d =", data_2d);
                PrintVector("pattern_data = ", pattern_data);
                PrintVV("pattern_data_2d =", pattern_data_2d);

                std::vector<float> similarity_result_temp;
                similarity_result_temp.resize(cols, 0);

                for (size_t i = 0; i < cols; i++)
                {
                    // normalize_vector(data_2d[i]);
                    // normalize_vector(pattern_data_2d[i]);
                    similarity_result_temp[i] = similarity_fun(data_2d[i], pattern_data_2d[i], distance_type);
                }
                PrintVector("similarity_result_temp = ", similarity_result_temp);

                similarity_result[pattern_index][0] = average_vector(similarity_result_temp);
            }
        }
        else
        {
            for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
            {
                std::vector<std::vector<float>> pattern_data_2d = ConvertVector1DTo2D(pattern_data[pattern_index], rows, cols, false);
                PrintVector("data = ", data);
                PrintVV("data_2d =", data_2d);
                PrintVector("pattern_data = ", pattern_data);
                PrintVV("pattern_data_2d =", pattern_data_2d);

                std::vector<float> similarity_result_temp;
                similarity_result_temp.resize(cols, 0);
                for (size_t shift_index = 0; shift_index < rows; shift_index++)
                {
                    for (size_t i = 0; i < cols; i++)
                    {
                        // normalize_vector(data_2d[i]);
                        // normalize_vector(pattern_data_2d[i]);
                        similarity_result_temp[i] = similarity_fun(data_2d[i], pattern_data_2d[i], distance_type);
                    }
                }
                PrintVector("similarity_result_temp = ", similarity_result_temp);
                similarity_result[pattern_index][0] = average_vector(similarity_result_temp);
            }
        }
    }
    else if (final_similarity_rank == 3)
    {
        int pattern_index = 0;
        std::vector<std::vector<std::vector<float>>> data_3d = ConvertVector1DTo3D(data, data_size[0], data_size[1], data_size[2], false);
        std::vector<std::vector<std::vector<float>>> pattern_data_3d = ConvertVector1DTo3D(pattern_data[pattern_index], each_pattern_data_size[0], each_pattern_data_size[1], each_pattern_data_size[2], false);
        std::vector<float> similarity_result_temp, similarity_result_temp2;
        size_t rows = data_size[0], cols = data_size[1], deps = data_size[2];
        similarity_result_temp.resize(cols, 0);
        similarity_result_temp2.resize(deps, 0);

        for (size_t j = 0; j < deps; j++)
        {
            for (size_t i = 0; i < cols; i++)
            {
                // normalize_vector(data_3d[j][i]);
                // normalize_vector(pattern_data_3d[j][i]);
                similarity_result_temp[i] = similarity_fun(data_3d[j][i], pattern_data_3d[j][i], distance_type);
            }
            similarity_result_temp2[j] = average_vector(similarity_result_temp);
        }
        similarity_result[pattern_index][0] = average_vector(similarity_result_temp2);
    }
    else
    {
        AU_EXIT("We haven't implement this 4D and beyond yet");
    }
    // similarity_result[pattern_index] = similarity_result_of_each_pattern;

    // return similarity_result;
}

template <class TT>
inline Stencil<std::vector<float>> tensor_search_udf(const Stencil<TT> &iStencil)
{
    std::vector<int> max_offset_upper, db_data_per_udf_size; // Size of input data for the whole chunk, zero based
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
    db_data_per_udf_size.resize(data_rank);
    for (int i = 0; i < data_rank; i++)
    {
        db_data_per_udf_size[i] = max_offset_upper[i] + 1;
    }
    // for (int pattern_index = 0; pattern_index < n_patterns; pattern_index++)
    // {
    //     similarity_vectors[pattern_index] = calculate_similarity(db_data_per_udf, db_data_per_udf_size, pattern_data[pattern_index], pattern_size);
    // }

    // calculate_similarity(std::vector<T1> &data, std::vector<int> &data_size, std::vector<std::vector<T2>> &pattern_data, std::vector<int> &each_pattern_data_size, int n_patterns, std::vector<std::vector<float>> &similarity_result)
    calculate_similarity(db_data_per_udf, db_data_per_udf_size, pattern_data, each_pattern_size, n_patterns, similarity_vectors);

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
    std::vector<int> pattern_chunk_size, pattern_overlap_size = {0, 0};
    std::vector<std::string> file_size_str_h5;

    std::string pattern_endpoint_id;
    bool is_pattern_input_single_file = false;

    if (isFile(extractDirFileName(pattern_dir)))
    {
        is_pattern_input_single_file = true;
        pattern_endpoint_id = db_data_file_type + ":" + pattern_dir;
    }
    else
    {
        pattern_endpoint_id = "EP_DIR:" + db_data_file_type + ":" + pattern_dir;
    }

    pattern_arrays = new AU::Array<float>(pattern_endpoint_id);

    if (!is_pattern_input_single_file)
    {
        pattern_arrays->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);
        pattern_arrays->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str_h5);
        std::cout << "file_size_str_h5 = " << file_size_str_h5[0] << "\n";
        String2Vector(file_size_str_h5[0], pattern_chunk_size);
    }
    else
    {
        std::vector<unsigned long long> array_size;
        pattern_arrays->GetArraySize(array_size);
        pattern_chunk_size.resize(array_size.size());
        for (int iii = 0; iii < array_size.size(); iii++)
        {
            pattern_chunk_size[iii] = array_size[iii];
        }
    }
    // assert(pattern_chunk_size.size() == data_rank);
    pattern_overlap_size.resize(pattern_chunk_size.size(), 0);
    pattern_arrays->SetChunkSize(pattern_chunk_size);
    pattern_arrays->SetOverlapSize(pattern_overlap_size);
    each_pattern_size = pattern_chunk_size;

    if (!ft_rank)
    {
        PrintVector("db-pattern: pattern_chunk_size = ", pattern_chunk_size);
        PrintVector("db-pattern: pattern_overlap_size = ", pattern_overlap_size);
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
    if (!is_pattern_input_single_file)
    {

        pattern_arrays->ControlEndpoint(DIR_N_FILES, file_size_str_h5);
        n_patterns = std::stoi(file_size_str_h5[0]);
    }
    else
    {
        n_patterns = 1;
    }

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

    // Special case, when user's data is stored as 2d but the similary happens on 1d
    if (ft_size == 1 && user_has_similarity_rank && pattern_chunk_size.size() == 2 && user_similarity_rank == 1)
    {
        n_patterns = pattern_chunk_size[0];
        int max_points_per_pattern = pattern_chunk_size[1];
        pattern_data = ConvertVector1DTo2D(pattern_data_per_file, n_patterns, max_points_per_pattern, true);
        each_pattern_size.resize(1);
        each_pattern_size[0] = pattern_chunk_size[1];
        // user_similarity_rank = data_rank;

        std::cout << "n_patterns = " << n_patterns << std::endl;
        std::cout << "each_pattern_size = " << each_pattern_size[0] << std::endl;
    }
}

int main(int argc, char *argv[])
{

    int copt;
    while ((copt = getopt(argc, argv, "d:q:r:m:s:k:h")) != -1)
        switch (copt)
        {
        case 'd':
            db_data_dir.assign(optarg);
            break;
        case 'q':
            pattern_dir.assign(optarg);
            break;
        case 'r':
            output_file_dir.assign(optarg);
            break;
        case 'm':
            distance_type.assign(optarg);
            break;
        case 's':
            is_shift = true;
            shift_size = atoi(optarg);
            std::cout << "Enable shift with shift_size = " << shift_size << "\n";
            break;
        case 'k':
            user_has_similarity_rank = true;
            user_similarity_rank = atoi(optarg);
            std::cout << "Enable the rank of similarity from user, user_similarity_rank = " << user_similarity_rank << "\n";
            break;
        case 'h':
            printf_help(argv[0]);
            exit(0);
        default:
            printf("Wrong option [%c] for %s \n", copt, argv[0]);
            printf_help(argv[0]);
            exit(-1);
            break;
        }

    // Init the MPICH, etc.
    AU_Init(argc, argv);

    if (isFile(extractDirFileName(db_data_dir)))
    {
        is_input_single_file = true;
        A_endpoint_id = db_data_file_type + ":" + db_data_dir;
    }
    else
    {
        A_endpoint_id = "EP_DIR:" + db_data_file_type + ":" + db_data_dir;
        is_input_single_file = false;
    }

    // A_endpoint_id = "EP_DIR:" + db_data_file_type + ":" + db_data_dir + ":" + db_data_h5_dataset;

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

    if (!is_input_single_file)
    {
        A->SkipFileTail();
        A->ExecuteUDFOnce();
        A->ControlEndpoint(DIR_SKIP_SIZE_CHECK, null_str);
    }

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

    if (!is_input_single_file)
    {
        std::vector<std::string> file_size_str;
        A->ControlEndpoint(DIR_GET_FILE_SIZE, file_size_str);
        String2Vector(file_size_str[0], chunk_size);
    }
    else
    {
        std::vector<unsigned long long> array_size;
        A->GetArraySize(array_size);
        chunk_size.resize(array_size.size());
        for (int iii = 0; iii < array_size.size(); iii++)
        {
            chunk_size[iii] = array_size[iii];
        }
    }

    data_rank = chunk_size.size();
    if (!user_has_similarity_rank)
    {
        user_similarity_rank = data_rank;
    }
    overlap_size.resize(data_rank, 0);
    A->SetChunkSize(chunk_size);
    A->SetOverlapSize(overlap_size);

    if (!ft_rank)
    {
        PrintVector("chunk_size = ", chunk_size);
        PrintVector("overlap_size = ", overlap_size);
    }

    AU::Array<float> *B;
    if (is_input_single_file)
    {
        // Store into a single file when the input file is a single file
        B = new AU::Array<float>(db_data_file_type + ":" + output_file_dir);
    }
    else
    {
        B = new AU::Array<float>("EP_DIR:" + db_data_file_type + ":" + output_file_dir);
        // Use the below rgx pattern to name the file
        std::vector<std::string> aug_output_replace_arg;
        aug_output_replace_arg.push_back(dir_output_match_rgx);
        aug_output_replace_arg.push_back(dir_output_replace_rgx);
        // B->ControlEndpoint(DIR_MERGE_INDEX, aug_merge_index);
        if (is_dir_output_match_replace_rgx)
            B->ControlEndpoint(DIR_OUPUT_REPLACE_RGX, aug_output_replace_arg);
    }

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

void printf_help(char *cmd)
{
    char *msg = (char *)"Usage: %s [OPTION]\n\
      	  -h help (--help)\n\
          -d string, [directory/file name]:[dataset name] of the database to search.\n\
          -q string, [directory/file name]:[dataset name]of query/pattern to search with.\n\
          -r string, [directory/file name]:[dataset name] to store result.\n\
          -m string, measurement of distance, dotproduct(default), euclidean, cosine, jaccard, angular,  \n\
          -s integer, enable shift with size on database (large tensor) to search, -1 means to shift by the size of pattern \n\
          Example: mpirun -n 1 %s -d db-data-2d:/testg/testd  -q db-pattern-2d:/testg/testd -r db-data-pattern-2d:/testg/testd\n";

    fprintf(stdout, msg, cmd, cmd);
}