#include <omp.h>
#include <zerg_file.h>
#include <zerg_template.h>
#include <zerg_time.h>
#include <regex>
#include <csv.h>

using namespace ztool;

struct ReduceMoment {
    void work() {
        if (threads == 0) threads = omp_get_max_threads();
        printf("thread num = %d\n", threads);
        omp_set_num_threads(threads);

        auto xs = ztool::path_wildcard(ztool::path_join(m_dir, "FeatureCorr.*.csv"));
        std::vector<std::pair<std::string, std::string>> todos;
        for (auto& item : xs) {
            int cob = std::stoi(item.first);
            if ((m_start_date > 0 && cob < m_start_date) || (m_end_date > 0 && cob > m_end_date)) {
                continue;
            } else {
                todos.emplace_back(item.first, item.second);
            }
        }
        std::sort(todos.begin(), todos.end(), [](auto& l, auto& r) { return l.first < r.first; });
        for (size_t i = 0; i < todos.size(); ++i) {
            auto& item = todos[i];
            printf("%s handle %s %zu/%zu\n", now_local_string().c_str(), item.second.c_str(), i, todos.size());
            std::cout << std::flush; // flush out
            reduce_single(item.first, item.second);
        }
        printf("finish\n");
    }

    std::string m_dir;
    int m_start_date{-1}, m_end_date{-1};
    int threads{0};
    std::vector<double> m_XTy;
    std::vector<std::vector<double>> m_XTX;

private:
    void reduce_single(const string& date_str, const string& path);
    size_t read_xx(const string& path, std::vector<std::vector<double>>& xx, std::vector<double>& mx);
};

static void help() {
    std::cout << "Program options:" << std::endl;
    std::cout << "  -h                                    list help" << std::endl;
    std::cout << "  -x arg (=)                          dir of x fst files" << std::endl;
    std::cout << "  -o arg (=)                          output dir / reduce input dir" << std::endl;
    printf("reduce_moment -x ~/junk/y_eval/data/ -o ~/junk/y_eval/output/\n");
}

int main(int argc, char** argv) {
    ReduceMoment dm;
    string config;
    int opt;
    while ((opt = getopt(argc, argv, "hx:o:")) != -1) {
        switch (opt) {
            case 'x':
                dm.m_dir = std::string(optarg);
                break;
            case 's':
                dm.m_start_date = std::stoi(optarg);
                break;
            case 'e':
                dm.m_end_date = std::stoi(optarg);
                break;
            case 'h':
            default:
                help();
                return 0;
        }
    }

    dm.work();
    return 0;
}

void ReduceMoment::reduce_single(const string& name, const string& path) {
    printf("reduce_single %s %s\n", name.c_str(), path.c_str());
    string xx_path = ztool::path_join(m_dir, "XX." + name + ".csv");
    if (!IsFileExisted(xx_path)) return;

    std::vector<std::vector<double>> xx;
    std::vector<double> mx;
    size_t valid_n = read_xx(xx_path, xx, mx);
    if (valid_n <= 0) return;

    std::vector<double> c_Y0_XY;
    std::vector<int> c_ukey, c_fid, c_Y0_N;

    io::CSVReader<4> infile(path);
    infile.read_header(io::ignore_extra_column, "ukey", "fid", "Y0_N", "Y0_XY");
    std::string fid;
    int ukey, Y0_N;
    double Y0_XY;
    while (infile.read_row(ukey, fid, Y0_N, Y0_XY)) {
        c_Y0_XY.push_back(Y0_XY);
        c_ukey.push_back(ukey);
        c_fid.push_back(std::stoi(fid.substr(1)));
        c_Y0_N.push_back(Y0_N);
    }

    printf("total row %zu\n", c_ukey.size());
    if (m_XTy.empty()) {
        m_XTy.resize(mx.size(), 0);
        m_XTX.resize(mx.size(), std::vector<double>(mx.size(), 0));
    }

    if (m_XTy.size() != mx.size()) {
        throw std::runtime_error("m_XTy size mismatch " + std::to_string(mx.size()));
    }
    if (m_XTX.size() != mx.size()) {
        throw std::runtime_error("m_XTX size mismatch " + std::to_string(mx.size()));
    }
    for (size_t i = 0; i < c_ukey.size(); ++i) {
        m_XTy[c_fid[i]] += c_Y0_XY[i];
    }
    for (size_t i = 0; i < mx.size(); ++i) {
        for (size_t j = 0; j < mx.size(); ++j) {
            m_XTX[i][j] += xx[i][j] * valid_n;
        }
    }
}

size_t ReduceMoment::read_xx(const string& path, std::vector<std::vector<double>>& xx, std::vector<double>& mx) {
    std::ifstream ifs(path);
    string line;
    ifs >> line;
    auto lets = ztool::split(line, ',');
    int N = (int)lets.size() - 1;
    mx = std::vector<double>(N, NAN);
    xx = std::vector<std::vector<double>>(N, std::vector<double>(N, NAN));
    size_t valid_n = std::stoul(lets[0]);
    for (int i = 0; i < N; ++i) {
        mx[i] = std::stod(lets[i + 1]);
        ifs >> line;
        auto lets1 = ztool::split(line, ',');
        for (int j = 0; j <= i; ++j) {
            xx[i][j] = std::stod(lets1[j]);
            if (i != j) xx[j][i] = xx[i][j];
        }
    }
    return valid_n;
}

