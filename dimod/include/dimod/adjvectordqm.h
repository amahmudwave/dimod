// Copyright 2020 D-Wave Systems Inc.
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#ifndef DIMOD_ADJVECTORDQM_H_
#define DIMOD_ADJVECTORDQM_H_

#include <stdio.h>
#include <algorithm>
#include <utility>
#include <vector>

#include "dimod/adjvectorbqm.h"
#include "dimod/utils.h"

namespace dimod {

template <class V, class B>
class AdjVectorDQM {
 public:
    using bias_type = B;
    using variable_type = V;
    using size_type = std::size_t;

    AdjVectorBQM bqm_;
    std::vector<variable_type> case_starts_;
    std::vector<std::vector<variable_type>> adj_;

    AdjVectorBQM() { case_starts_.push_back(0); }

    variable_type num_variables() { return adj_.size(); }

    size_type num_variaables_interactions() {
        size_type num = 0;
        for (auto v = 0, vend = this->num_variables(); v < vend; v++) {
            num += adj_[v].size();
        }
        return (num / 2);
    }

    variable_type num_cases(variable_type v = -1) {
        assert(v < this->num_variables());
        if (v < 0) {
            return bqm_.num_variables();
        } else {
            return (case_starts_[v + 1] - case_starts_[v]);
        }
    }

    size_type num_case_interactions() { return bqm_.num_interactions(); }

    variable_type add_variable(size_t num_cases) {
        assert(num_cases > 0);
        auto v = adj_.size();
        adj_.resize(v + 1);
        for (auto n = 0; n < num_cases; n++) {
            bqm_.add_variable();
        }
        case_starts_.push_back(bqm_.num_variables());
        return v;
    }

    void get_linear(variable_type v, bias_type* biases) {
        assert(v >= 0 && v < this->num_variables());
        for (auto case_v = 0, num_cases_v = this->num_cases(v);
             case_v < num_cases_v; case_v++) {
            biases[case_v] = bqm_.get_linear(case_starts_[v] + case_v);
        }
    }

    bias_type get_linear_case(variable_type v, variable_type case_v) {
        assert(v >= 0 && v < this->num_variables());
        assert(case_v >= 0 && case_v < num_cases(v));
        return bqm_.get_linear(case_starts_[v] + case_v);
    }

    void set_linear(variable_type v, bias_type* p_biases) {
        for (auto case_v = 0, num_cases_v = this->num_cases(v);
             case_v < num_cases_v; case_v++) {
            bqm_.set_linear(case_starts_[v] + case_v, p_biases[case_v]);
        }
    }

    void set_linear_case(variable_type v, variable_type case_v, bias_type b) {
        assert(case_v >= 0 && case_v < this->num_cases(v));
        bqm_.set_linear(case_starts_[v] + case_v, b);
    }

    // Returns false if there is no interaction among the variables.
    bool get_quadratic(variable_type u, variable_type v,
                       bias_type* quadratic_biases) {
        assert(u >= 0 && u < this->num_variables());
        assert(v >= 0 && v < this->num_variables());
        auto it = std::lower_bound(adj_[u].begin(), adj_[u].end(), v);
        if (it == adj_[u].end() || *it != v) {
            return false;
        }
        auto num_cases_u = num_cases(u);
        auto num_cases_v = num_cases(v);
        for (auto case_u = 0; case_u < num_cases_u; case_u++) {
            auto span = bqm_.neighborhood(case_starts_[u] + case_u,
                                          case_starts_[v]);
            while (span.first != span.second &&
                   *(span.first) < case_starts_[v + 1]) {
                case_v = *(span.first) - case_starts_[v];
                quadratic_biases[case_u * num_cases_v + case_v] =
                        *(span.first).second;
                span.first++;
            }
        }
        return true;
    }

    bias_type get_quadratic_case(variable_type u, variable_type case_u,
                                 variable_type v, variable_type case_v) {
        assert(u >= 0 && u < this->num_variables());
        assert(case_u >= 0 && case_u < num_cases(v));
        assert(v >= 0 && v < this->num_variables());
        assert(case_v >= 0 && case_v < num_cases(v));
        // should add assert for u != v ?
        auto cu = case_starts_[u] + case_u;
        auto cv = case_starts_[v] + case_v;
        return bqm_.get_quadratic(cu, cv).first;
    }

    bool set_quadratic(variable_type u, variable_type v, bias_type* p_biases) {
        assert(u >= 0 && u < this->num_variables());
        assert(v >= 0 && v < this->num_variables());
        assert(u != v);
        num_cases_u = num_cases(u);
        num_cases_v = num_cases(v);
        auto num_cases_u = num_cases(u);
        auto num_cases_v = num_cases(v);
        for (auto case_u = 0; case_u < num_cases_u; case_u++) {
            cu = case_starts_[u] + case_u;
            for (auto case_v = 0; case_v < num_cases_v; case_v++) {
                cv = case_starts_[v] + case_v;
                auto bias = p_biases[cu * num_cases_v + case_v];
                bqm_.set_quadratic(cu, cv, bias);
            }
        }
        auto low = std::lower_bound(adj_[u].begin(), adj_[u].end(), v);
        if (low == adj_[u].end() || *low != v) {
            adj_[u].insert(low, v);
            adj_[v].insert(std::lower_bound(adj_[v].begin(), adj_[v].end(), u),
                           u);
        }
        return true;
    }

    // Check if boolean type is still okay
    bool set_quadratic_case(variable_type u, variable_type case_u,
                            variable_type v, variable_type case_v,
                            bias_type bias) {
        assert(u >= 0 && u < this->num_variables());
        assert(case_u >= 0 && case_u < num_cases(v));
        assert(v >= 0 && v < this->num_variables());
        assert(case_v >= 0 && case_v < num_cases(v));
        auto cu = case_starts_[u] + case_u;
        auto cv = case_starts_[v] + case_v;
        bqm_.set_quadratic(cu, cv, bias);
        auto low = std::lower_bound(adj_[u].begin(), adj_[u].end(), v);
        if (low == adj_[u].end() || *low != v) {
            adj_[u].insert(low, v);
            adj_[v].insert(std::lower_bound(adj_[v].begin(), adj_[v].end(), u),
                           u);
        }
        return true;
    }

    void energies(variable_type* p_samples, int num_samples,
                  variable_type num_variables, bias_type* p_energies) {
        assert(num_variables == this->num_variables());
        memset(p_energies, 0, num_samples * sizeof(bias_type));
#pragma omp parallel for
        for (auto si = 0; si < num_samples; si++) {
            variable_type* p_curr_sample = samples + (si * num_variables);
            bias_type* p_curr_energy = p_energies + si;
            for (auto u = 0; u < num_variables; u++) {
                auto case_u = p_curr_sample[u];
                assert(case_u < num_cases(u));
                auto cu = case_starts_[u] + case_u;
                *p_curr_energy += bqm_.get_linear(cu);
                for (auto vi = 0; vi < adj_[u].size(); vi++) {
                    auto v = adj_[u][vi];
                    // We only care about lower triangle.
                    if (v > u) {
                        break;
                    }
                    auto case_v = p_cur_sample[v];
                    auto cv = case_starts_[v] + case_v;
                    auto out = bqm_.get_quadratic(cu, cv);
                    if (out.second) {
                        p_energies[si] += out.first;
                    }
                }
            }
        }
    }
}
}  // namespace dimod

#endif  // DIMOD_ADJVECTORDQM_H_