#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iostream>
#include <libconfig.h++>
#include <string>
#include <unistd.h>
#include <vector>

bool check_unique_mutation(const int ind, 
                           const std::vector<std::vector<int> > &hap) {
    int ret = 1;
    int j = 0;
    while (j < (int)hap.size() && ret == 1) {
        if (j != ind && std::binary_search(hap[j].begin(), hap[j].end(), 
                                            hap[ind].back())) {
            ret = 0;
        }
        j++;
    }
    return ret;
}

double coal_event(const int island, const int ind, 
                  std::vector<std::vector<int> > &count, const double K, 
                  const double proposal, const std::vector<double> &n, 
                  const double r) {
    double log_w = log(n[island]) + log((double)(count[island][ind] - 1)) 
        - log(r) - log(proposal);
    if (island == 1) {
        log_w += log(K);
    }
    count[island][ind]--;
    return log_w;
}

double mut_event(const int island, const int model, const int ind,
                 std::vector<std::vector<int> > &hap,
                 std::vector<std::vector<int> > &count, const double theta,
                 const double theta_bank, const double proposal, 
                 const double r) {
    hap[ind].pop_back();
    double log_w = -log(proposal) - log(r);
    if (island == 1 && model == 2) {
        log_w += log(theta_bank);
    } else {
        log_w += log(theta);
    }
    int hap_index = 0;
    int check = 0;
    while (check == 0 && hap_index < (int)hap.size()) {
        if (ind != hap_index && hap[ind] == hap[hap_index]) {
            check = 1;
        } else {
            hap_index++;
        }
    }
    if (check == 1) {
        log_w += log((double)(count[island][hap_index] + 1));
        count[island][hap_index]++;
        hap.erase(hap.begin() + ind);
        count[0].erase(count[0].begin() + ind);
        count[1].erase(count[1].begin() + ind);
    }
    return log_w;
}

double mig_event(const int island, const int ind,
                 std::vector<std::vector<int> > &count, const double c, 
                 const double K, const std::vector<double> &n, 
                 const double proposal, const double r) {
    double log_w = log(c) + log(n[island]) 
        + log((double)(count[(island + 1) % 2][ind] + 1)) 
        - log(n[(island + 1) % 2] + 1.0) - log(proposal) - log(r);
    if (island == 1) {
        log_w += log(K);
    }
    count[island][ind]--;
    count[(island + 1) % 2][ind]++;
    return log_w;
}

void assign_probabilities_kingman(double &coal, double &mut, 
                                  const int island, const int ind,
                                  const std::vector<std::vector<int> > &count,
                                  const std::vector<std::vector<int> > &hap) {
    if (count[island][ind] == 1 && check_unique_mutation(ind, hap) == 1) {
        mut = 1.0;
    } else if (count[island][ind] > 1) {
        coal = 1.0;
    }
    return;
}

void assign_probabilities_strong(double &coal, double &mig, double &mut, 
                                 const int island, const int ind,
                                 const std::vector<std::vector<int> > &count,
                                 const std::vector<std::vector<int> > &hap,
                                 const double theta, const double theta_bank,
                                 const double c, const double K) {
    coal = fmax((double)(count[0][ind] - 1), 0.0);
    mig = c;
    if (count[0][ind] + count[1][ind] == 1 
        && check_unique_mutation(ind, hap) == 1) {
        mut = theta;
    }
    if (island == 1) {
        coal = 0.0;
        mig *= K;
        if (count[0][ind] + count[1][ind] == 1 
            && check_unique_mutation(ind, hap) == 1) {
            mut = theta_bank;
        }
    }
    return;
}

void assign_probabilities_structured(double &coal, double &mig, double &mut, 
                                     const int island, const int ind,
                                     const std::vector<std::vector<int> > &count,
                                     const std::vector<std::vector<int> > &hap,
                                     const std::vector<double> &n, 
                                     const double theta, const double c, 
                                     const double K) {
    double r_0 = fmax(n[0] - 1.0, 0.0) + c + theta;
    double r_1 = fmax(n[1] - 1.0, 0.0) * K + K * c + theta;
    double pi_0 = (r_1 * fmax((double)(count[0][ind] - 1), 0.0) + c * K 
        * fmax((double)(count[1][ind] - 1), 0.0)) / (r_0 * r_1 - K * c * c);
    double pi_1 = (K * c * fmax((double)(count[0][ind] - 1), 0.0) + r_1 * K
        * fmax((double)(count[1][ind] - 1), 0.0)) / (r_0 * r_1 - K * c * c);
    if (n[0] + n[1] == 2) {
        pi_0 = 1.0;
        pi_1 = 1.0;
    }
    if (island == 0) {
        if (count[0][ind] + count[1][ind] == 1 
        && check_unique_mutation(ind, hap) == 1) {
            mig = c;
            mut = theta;
        } else if (count[0][ind] > 1) {
            coal = (double)(count[0][ind] - 1);
            mig = c * pi_1;
        } else {
            mig = 1.0;
        }
    } else {
        if (count[0][ind] + count[1][ind] == 1 
        && check_unique_mutation(ind, hap) == 1) {
            mig = c * K;
            mut = theta;
        } else if (count[1][ind] > 1) {
            coal = K * (double)(count[1][ind] - 1);
            mig = K * c * pi_0;
        } else {
            mig = 1.0;
        }
    }
    return;
}

double step(std::vector<std::vector<int> > &hap, 
            std::vector<std::vector<int> > &count, const double theta, 
            const double theta_bank, const double c, const double K, 
            gsl_rng *gen, int &total_count, const int model) {
    std::vector<double> n(2, 0.0);
    // For Kingman when not every lineage can necessarily be hit by an event
    double n_prop = 0.0; 
    for (unsigned int i = 0; i < hap.size(); i++) {
        n[0] += (double)count[0][i];
        n[1] += (double)count[1][i];
        if (model == 0) {
            if (count[0][i] == 1 
                && check_unique_mutation(i, hap) == 1) {
                n_prop += 1.0;
            } else if (count[0][i] > 1) {
                n_prop += (double)count[0][i];
            }
        }
    }
    std::vector<double> r = n;
    switch (model) {
        case 0: r[0] *= r[0] - 1.0 + theta;
                break;
        case 1: r[0] *= fmax(r[0] - 1.0, 0.0) + c + theta;
                r[1] *= fmax(r[1] - 1.0, 0.0) * K + K * c + theta;
                break;
        case 2: r[0] *= fmax(r[0] - 1.0, 0.0) + c + theta;
                r[1] *= K * c + theta_bank;
                break;
    }
    int island = 0;
    if (gsl_rng_uniform(gen) < r[1] / (r[0] + r[1])) {
        island = 1;
    }
    double coin = gsl_rng_uniform_pos(gen);
    int ind = -1;
    double ub = 0.0;
    while (coin > ub) {
        ind++;
        if (model == 0) {
            if (count[0][ind] == 1 && check_unique_mutation(ind, hap) == 1) {
                ub += 1.0 / n_prop;
            } else if (count[0][ind] > 1) {
                ub += (double)count[0][ind] / n_prop;
            }
        } else {
            ub += (double)count[island][ind] / n[island];
        }
    }
    double coal = 0.0, mig = 0.0, mut = 0.0;
    switch (model) {
        case 0: assign_probabilities_kingman(coal, mut, island, ind, count, 
                                             hap);
                break;
        case 1: assign_probabilities_structured(coal, mig, mut, island, ind, 
                                                count, hap, n, theta, c, K);
                break;
        case 2: assign_probabilities_strong(coal, mig, mut, island, ind, count, 
                                            hap, theta, theta_bank, c, K);
                break;
    }
    coin = gsl_rng_uniform(gen);
    double total = coal + mut + mig;
    double log_w = 0.0;
    if (coin < coal / total) {
        if (model == 0) {
            log_w = coal_event(0, ind, count, K, 
                               (double)count[0][ind] / n_prop, n, r[0]);
        } else {
            log_w = coal_event(island, ind, count, K, coal * r[island] 
                * (double)count[island][ind] / (total * n[island] 
                * (r[0] + r[1])), n, r[0] + r[1]);
        }
        total_count--;
    } else if (coin < (coal + mig) / total) {
        log_w = mig_event(island, ind, count, c, K, n, mig * r[island] 
            * (double)count[island][ind] / (total * n[island] * (r[0] + r[1])), 
            r[0] + r[1]);
    } else {
        if (model == 0) {
            log_w = mut_event(0, 0, ind, hap, count, theta, theta_bank, 
                              (double)count[0][ind] / n_prop, r[0]);
        } else {
            log_w = mut_event(island, model, ind, hap, count, theta, theta_bank, 
                mut * r[island] * (double)count[island][ind] / (total 
                * n[island] * (r[0] + r[1])), r[0] + r[1]);
        }
    }
    return log_w;
}

double smc(const int model, const double theta, const double theta_bank, 
           const double K, const double c, const int particle_number, 
           const std::vector<std::vector<int> > &count,
           const std::vector<std::vector<int> > &hap, gsl_rng *gen) {
    std::vector<std::vector<std::vector<int> > > N_hap(particle_number, hap);
    std::vector<std::vector<std::vector<int> > > N_count(particle_number, 
                                                         count);
    std::vector<double> log_w(particle_number, 0.0);
    std::vector<int> anc(particle_number);
    double log_sum_w = log((double)particle_number);
    double log_sum_w_squared = log((double)particle_number);
    int limit = 0;
    for (unsigned int i = 0; i < count[0].size(); i++) {
        limit += count[0][i] + count[1][i];
    }
    limit--;
    double coin, rhs, max_log_w;
    int anc_ind, total_count;
    while (limit > 0) {
        if (2.0 * log_sum_w - log_sum_w_squared < log((double)particle_number 
            / 2.0)) {
            anc_ind = 0;
            rhs = log_w[0] - log_sum_w;
            for (int i = 0; i < particle_number; i++) {
                coin = log(((double)i + gsl_rng_uniform_pos(gen)) 
                    / (double)particle_number);
                while (coin > rhs) {
                    anc_ind++;
                    rhs = log(exp(rhs) + exp(log_w[anc_ind] - log_sum_w));
                }
                anc[i] = anc_ind;
            }
            int i_tmp = -1;
            for (int i = 0; i < particle_number; i++) {
                if (anc[i] != i && anc[anc[i]] != anc[i]) {
                    i_tmp = anc[anc[i]];
                    anc[anc[i]] = anc[i];
                    anc[i] = i_tmp;
                    i--;
                }
            }
            for (int i = 0; i < particle_number; i++) {
                if (anc[i] != i) {
                    N_hap[i] = N_hap[anc[i]];
                    N_count[i] = N_count[anc[i]];
                }
                log_w[i] = log_sum_w - log((double)particle_number);
            }
        }
        for (int i = 0; i < particle_number; i++) {
            total_count = limit + 1;
            while (total_count > limit) {
                log_w[i] += step(N_hap[i], N_count[i], theta, theta_bank, c, K,
                                 gen, total_count, model);
            }
        }
        limit--;
        max_log_w = log_w[0];
        for (int i = 1; i < particle_number; i++) {
            max_log_w = fmax(max_log_w, log_w[i]);
        }
        log_sum_w = exp(log_w[0] - max_log_w);
        log_sum_w_squared = exp(2.0 * (log_w[0] - max_log_w));
        for (int i = 1; i < particle_number; i++) {
            log_sum_w += exp(log_w[i] - max_log_w);
            log_sum_w_squared += exp(2.0 * (log_w[i] - max_log_w));
        }
        log_sum_w = log(log_sum_w) + max_log_w;
        log_sum_w_squared = log(log_sum_w_squared) + 2.0 * max_log_w;
    }
    return log_sum_w - log((double)particle_number);
}
