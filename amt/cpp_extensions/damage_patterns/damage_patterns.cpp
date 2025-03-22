#include <Python.h>
#include <numpy/arrayobject.h>
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>
#include "damage_patterns.h"

struct DamageStats {
    std::vector<double> c_to_t_5prime;
    std::vector<double> g_to_a_3prime;
    std::vector<int> fragment_lengths;
    
    DamageStats(int window_size, int max_length) {
        c_to_t_5prime.resize(window_size, 0.0);
        g_to_a_3prime.resize(window_size, 0.0);
        fragment_lengths.resize(max_length + 1, 0);
    }
};

class DamageCounter {
private:
    std::mutex mtx;
    int window_size;
    int max_length;
    
    // Counts for different mismatches at each position
    std::vector<int> c_t_5prime_counts;
    std::vector<int> g_a_3prime_counts;
    std::vector<int> c_5prime_total;
    std::vector<int> g_3prime_total;
    std::vector<int> fragment_length_counts;

public:
    DamageCounter(int ws, int ml) : window_size(ws), max_length(ml) {
        c_t_5prime_counts.resize(window_size, 0);
        g_a_3prime_counts.resize(window_size, 0);
        c_5prime_total.resize(window_size, 0);
        g_3prime_total.resize(window_size, 0);
        fragment_length_counts.resize(max_length + 1, 0);
    }
    
    void addRead(int length, const std::vector<int>& ct_5p, const std::vector<int>& ga_3p,
                 const std::vector<int>& c_5p, const std::vector<int>& g_3p) {
        std::lock_guard<std::mutex> lock(mtx);
        
        if (length >= 0 && length <= max_length) {
            fragment_length_counts[length]++;
        }
        
        for (int i = 0; i < window_size; i++) {
            c_t_5prime_counts[i] += ct_5p[i];
            g_a_3prime_counts[i] += ga_3p[i];
            c_5prime_total[i] += c_5p[i];
            g_3prime_total[i] += g_3p[i];
        }
    }
    
    void merge(const DamageCounter& other) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (int i = 0; i <= max_length; i++) {
            fragment_length_counts[i] += other.fragment_length_counts[i];
        }
        
        for (int i = 0; i < window_size; i++) {
            c_t_5prime_counts[i] += other.c_t_5prime_counts[i];
            g_a_3prime_counts[i] += other.g_a_3prime_counts[i];
            c_5prime_total[i] += other.c_5prime_total[i];
            g_3prime_total[i] += other.g_3prime_total[i];
        }
    }
    
    DamageStats getStats() const {
        DamageStats stats(window_size, max_length);
        
        for (int i = 0; i < window_size; i++) {
            if (c_5prime_total[i] > 0) {
                stats.c_to_t_5prime[i] = static_cast<double>(c_t_5prime_counts[i]) / c_5prime_total[i];
            }
        }
        
        for (int i = 0; i < window_size; i++) {
            if (g_3prime_total[i] > 0) {
                stats.g_to_a_3prime[i] = static_cast<double>(g_a_3prime_counts[i]) / g_3prime_total[i];
            }
        }
        
        // Copy fragment length distribution
        stats.fragment_lengths = fragment_length_counts;
        
        return stats;
    }
};

void processRead(bam1_t* read, const char* ref_seq, int ref_start, int window_size,
                std::vector<int>& ct_5p, std::vector<int>& ga_3p,
                std::vector<int>& c_5p, std::vector<int>& g_3p) {
    
    int read_len = read->core.l_qseq;
    uint8_t* qseq = bam_get_seq(read);
    uint32_t* cigar = bam_get_cigar(read);
    uint32_t n_cigar = read->core.n_cigar;
    
    int read_pos = 0;
    int ref_pos = 0;
    
    for (uint32_t i = 0; i < n_cigar; i++) {
        uint32_t op = cigar[i] & BAM_CIGAR_MASK;
        uint32_t len = cigar[i] >> BAM_CIGAR_SHIFT;
        
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (uint32_t j = 0; j < len; j++) {
                char ref_base = ref_seq[ref_pos + ref_start];
                char read_base = seq_nt16_str[bam_seqi(qseq, read_pos)];
                
                // Check 5' end for C→T
                if (read_pos < window_size) {
                    if (ref_base == 'C') {
                        c_5p[read_pos]++;
                        if (read_base == 'T') {
                            ct_5p[read_pos]++;
                        }
                    }
                }
                
                // Check 3' end for G→A
                if (read_pos >= read_len - window_size) {
                    int pos_from_end = read_len - read_pos - 1;
                    if (ref_base == 'G') {
                        g_3p[pos_from_end]++;
                        if (read_base == 'A') {
                            ga_3p[pos_from_end]++;
                        }
                    }
                }
                
                read_pos++;
                ref_pos++;
            }
        } else if (op == BAM_CINS) {
            read_pos += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            ref_pos += len;
        } else if (op == BAM_CSOFT_CLIP) {
            read_pos += len;
        }
    }
}

// Process a chunk of BAM file in a thread
void processBamChunk(const char* bam_path, const char* ref_path, 
                     const std::string& region, int window_size, 
                     int min_length, int max_length, DamageCounter& counter) {
    
    samFile* bam = sam_open(bam_path, "r");
    if (!bam) return;
    
    bam_hdr_t* header = sam_hdr_read(bam);
    if (!header) {
        sam_close(bam);
        return;
    }
    
    hts_idx_t* idx = sam_index_load(bam, bam_path);
    if (!idx) {
        bam_hdr_destroy(header);
        sam_close(bam);
        return;
    }
    
    // Open reference
    faidx_t* fai = fai_load(ref_path);
    if (!fai) {
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(bam);
        return;
    }
    
    // Set up an iterator for the region
    hts_itr_t* iter = NULL;
    if (region.empty()) {
        iter = sam_itr_queryi(idx, HTS_IDX_START, 0, 0);
    } else {
        iter = sam_itr_querys(idx, header, region.c_str());
    }
    
    if (!iter) {
        fai_destroy(fai);
        hts_idx_destroy(idx);
        bam_hdr_destroy(header);
        sam_close(bam);
        return;
    }
    
    bam1_t* read = bam_init1();
    
    std::vector<int> ct_5p(window_size, 0);
    std::vector<int> ga_3p(window_size, 0);
    std::vector<int> c_5p(window_size, 0);
    std::vector<int> g_3p(window_size, 0);
    
    while (sam_itr_next(bam, iter, read) >= 0) {
        if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FDUP)) {
            continue;
        }
        
        int qlen = read->core.l_qseq;
        if (qlen < min_length || qlen > max_length) {
            continue;
        }
        
        int tid = read->core.tid;
        const char* ref_name = header->target_name[tid];
        int ref_start = read->core.pos;
        int ref_end = bam_endpos(read);
        
        int seq_len;
        char* ref_seq = faidx_fetch_seq(fai, ref_name, ref_start, ref_end - 1, &seq_len);
        if (!ref_seq || seq_len <= 0) {
            if (ref_seq) free(ref_seq);
            continue;
        }
        
        std::fill(ct_5p.begin(), ct_5p.end(), 0);
        std::fill(ga_3p.begin(), ga_3p.end(), 0);
        std::fill(c_5p.begin(), c_5p.end(), 0);
        std::fill(g_3p.begin(), g_3p.end(), 0);
        
        processRead(read, ref_seq, 0, window_size, ct_5p, ga_3p, c_5p, g_3p);
        
        counter.addRead(qlen, ct_5p, ga_3p, c_5p, g_3p);
        
        free(ref_seq);
    }
    
    bam_destroy1(read);
    hts_itr_destroy(iter);
    fai_destroy(fai);
    hts_idx_destroy(idx);
    bam_hdr_destroy(header);
    sam_close(bam);
}

// Split BAM file into regions for parallel processing
std::vector<std::string> splitBamRegions(const char* bam_path, int num_threads) {
    std::vector<std::string> regions;
    
    samFile* bam = sam_open(bam_path, "r");
    if (!bam) return regions;
    
    bam_hdr_t* header = sam_hdr_read(bam);
    if (!header) {
        sam_close(bam);
        return regions;
    }
    
    // Regions based on chromosomes
    for (int i = 0; i < header->n_targets; i++) {
        regions.push_back(std::string(header->target_name[i]));
    }
    
    bam_hdr_destroy(header);
    sam_close(bam);
    
    return regions;
}

// Python wrapper function
static PyObject* analyze_damage_wrapper(PyObject* self, PyObject* args) {
    const char* bam_path;
    const char* ref_path;
    int min_length, max_length, window_size, threads;
    
    if (!PyArg_ParseTuple(args, "ssiiii", &bam_path, &ref_path, 
                         &min_length, &max_length, &window_size, &threads)) {
        return NULL;
    }
    
    // Get regions for parallel processing
    std::vector<std::string> regions = splitBamRegions(bam_path, threads);
    if (regions.empty()) {
        PyErr_SetString(PyExc_RuntimeError, "Failed to split BAM into regions");
        return NULL;
    }
    
    DamageCounter master_counter(window_size, max_length);
    std::vector<DamageCounter> thread_counters;
    std::vector<std::thread> workers;
    
    // Determine actual thread count
    int actual_threads = std::min(threads, static_cast<int>(regions.size()));
    if (actual_threads <= 0) actual_threads = 1;
    for (int i = 0; i < actual_threads; i++) {
        thread_counters.emplace_back(window_size, max_length);
    }
    for (int t = 0; t < actual_threads; t++) {
        workers.emplace_back([&, t]() {
            for (size_t i = t; i < regions.size(); i += actual_threads) {
                processBamChunk(bam_path, ref_path, regions[i], 
                               window_size, min_length, max_length, thread_counters[t]);
            }
        });
    }
    
    for (auto& worker : workers) {
        worker.join();
    }
    
    for (const auto& counter : thread_counters) {
        master_counter.merge(counter);
    }
    
    DamageStats stats = master_counter.getStats();
    
    PyObject* result_dict = PyDict_New();
    
    import_array();
    
    npy_intp fragment_dims[1] = {static_cast<npy_intp>(max_length + 1)};
    PyObject* fragment_array = PyArray_SimpleNew(1, fragment_dims, NPY_INT);
    std::memcpy(PyArray_DATA((PyArrayObject*)fragment_array), 
                stats.fragment_lengths.data(), 
                stats.fragment_lengths.size() * sizeof(int));
    PyDict_SetItemString(result_dict, "fragment_length", fragment_array);
    
    npy_intp window_dims[1] = {static_cast<npy_intp>(window_size)};
    PyObject* ct_5prime = PyArray_SimpleNew(1, window_dims, NPY_DOUBLE);
    std::memcpy(PyArray_DATA((PyArrayObject*)ct_5prime), 
                stats.c_to_t_5prime.data(), 
                window_size * sizeof(double));
    
    PyObject* ga_3prime = PyArray_SimpleNew(1, window_dims, NPY_DOUBLE);
    std::memcpy(PyArray_DATA((PyArrayObject*)ga_3prime), 
                stats.g_to_a_3prime.data(), 
                window_size * sizeof(double));
    
    PyObject* five_prime_dict = PyDict_New();
    PyDict_SetItemString(five_prime_dict, "C→T", ct_5prime);
    
    PyObject* three_prime_dict = PyDict_New();
    PyDict_SetItemString(three_prime_dict, "G→A", ga_3prime);
    
    PyDict_SetItemString(result_dict, "5_prime", five_prime_dict);
    PyDict_SetItemString(result_dict, "3_prime", three_prime_dict);
    
    return result_dict;
}

static PyMethodDef DamageMethods[] = {
    {"analyze_damage", analyze_damage_wrapper, METH_VARARGS,
     "Analyze ancient DNA damage patterns"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef damagemodule = {
    PyModuleDef_HEAD_INIT,
    "damage_patterns",
    "C++ extension for ancient DNA damage pattern analysis",
    -1,
    DamageMethods
};

PyMODINIT_FUNC PyInit_damage_patterns(void) {
    import_array();
    return PyModule_Create(&damagemodule);
}
