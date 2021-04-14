#ifndef ORNATE_MATH_GRAPH_H
#define ORNATE_MATH_GRAPH_H

#include <algorithm>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct GraphSplitter {
    GraphSplitter(int max_node_cnt, int max_group_node_cnt) {
        m_max_node_cnt = max_node_cnt;
        m_max_group_nodes = max_group_node_cnt;
        m_node2cnt.resize(m_max_node_cnt, 0);
        m_node2constrains.resize(m_max_node_cnt, std::vector<int>());
    }

    /**
     * @param constrain 包含 node index, 应该是个set, node不重复
     */
    void add_constrain(const std::vector<int>& constrain) {
        m_constrains.push_back(constrain);
        int ci = (int)m_constrains.size() - 1;
        for (int ni : constrain) {
            m_node2cnt[ni]++;
            m_node2constrains[ni].push_back(ci);
        }
    }
    bool work() {
        total_constrains = (int)m_constrains.size();
        m_result_constrains.clear();
        m_result_nodes.clear();
        filled.resize(total_constrains);
        std::fill(filled.begin(), filled.end(), false);

        while (init_round()) {
            while (!m_node_q.empty()) {
                int ni = m_node_q.top();
                m_node_q.pop();

                const auto& constrains_ = m_node2constrains[ni];
                for (int ci : constrains_) {
                    if (try_add_constrain(ci)) break;
                }
            }

            // re-insure no missing for curr round
            for (int ci = 0; ci < total_constrains; ++ci) {
                try_add_constrain(ci);
            }

            if (!check_curr_round_finish()) return false;
        }
        return true;
    }

    std::vector<int> get_nodes(int idx) {
        const std::unordered_set<int>& s = m_result_nodes[idx];
        return std::vector<int>(s.begin(), s.end());
    }

private:
    bool init_round() {
        m_curr_constrains.clear();
        m_curr_nodes.clear();
        m_node_q = std::priority_queue<int>();
        std::unordered_map<int, int> node2cnt;
        for (int i = 0; i < total_constrains; ++i) {
            if (filled[i]) continue;

            for (int di : m_constrains[i]) {
                node2cnt[di] = m_node2cnt[di];
            }
        }

        if (node2cnt.empty()) return false;
        std::vector<std::pair<int, int>> sort_;
        std::transform(node2cnt.begin(), node2cnt.end(), std::back_inserter(sort_), [](auto& item) { return item; });
        std::sort(sort_.begin(), sort_.end(), [](auto& l, auto& r) { return l.second < r.second; });
        for (const auto& item : sort_) m_node_q.push(item.first);
        return true;
    }

    bool try_add_constrain(int ci) {
        if (filled[ci]) return false;

        int need_add_node_cnt = 0;
        for (int ni : m_constrains[ci]) {
            if (m_curr_nodes.find(ni) == m_curr_nodes.end()) {
                ++need_add_node_cnt;
            }
        }

        if (need_add_node_cnt <= m_max_group_nodes - (int)m_curr_nodes.size()) {
            for (int ni : m_constrains[ci]) {
                m_curr_nodes.insert(ni);
            }
            m_curr_constrains.push_back(ci);
            filled[ci] = true;
            // printf("add ci=%d, nodes=%s\n", ci, head(m_constrains[ci], 0).c_str());
            return true;
        }
        return false;
    }
    bool check_curr_round_finish() {
        if (m_curr_constrains.empty()) return false;
        m_result_constrains.push_back(m_curr_constrains);
        m_result_nodes.push_back(m_curr_nodes);
        // printf("archive curr %s\n", head(m_curr_constrains, 0).c_str());
        return true;
    }

public:
    int total_constrains{0};
    int m_max_node_cnt{0};
    int m_max_group_nodes{0};
    std::vector<bool> filled;
    std::vector<int> m_curr_constrains;
    std::unordered_set<int> m_curr_nodes;
    std::priority_queue<int> m_node_q;
    std::vector<int> m_node2cnt;
    std::vector<std::vector<int>> m_constrains;  // constrain means one entry must be together
    std::vector<std::vector<int>> m_node2constrains;
    std::vector<std::vector<int>> m_result_constrains;
    std::vector<std::unordered_set<int>> m_result_nodes;
};

#endif
