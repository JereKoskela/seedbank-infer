#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_gamma.h>
#include <iomanip>
#include <iostream>
#include <libconfig.h++>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>
#include "smc.hh"

int main (int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Call " << argv[0] << " <config file>" << std::endl;
        return 1;
    }
    gsl_rng *gen = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(gen, time(NULL) * getpid());
    libconfig::Config cfg;
    cfg.readFile(argv[1]);
    std::vector<int> models;
    std::vector<int> mutation_flags;
    for (int i = 0; i < cfg.getRoot()["models"].getLength(); i++) {
        models.push_back(cfg.getRoot()["models"][i]);
        mutation_flags.push_back(cfg.getRoot()["mutation_flags"][i]);
    }
    std::vector<int> particle_numbers;
    for (int i = 0; i < cfg.getRoot()["particle_numbers"].getLength(); i++) {
        particle_numbers.push_back(cfg.getRoot()["particle_numbers"][i]);
    }
    int run_length;
    cfg.lookupValue("run_length", run_length);
    double sd;
    cfg.lookupValue("st_dev", sd);
    std::string data_file_str;
    cfg.lookupValue("data_file", data_file_str);
    std::ifstream data_file;
    data_file.open(data_file_str);
    if (!data_file.good()) {
        std::cerr << "invalid data file path." << std::endl;
        return 1;
    }
    std::string line, token;
    std::stringstream iss;
    std::vector<int> tmp;
    std::vector<std::vector<int> > hap;
    std::vector<std::vector<int> > count(2);
    getline(data_file, line);
    while (!data_file.eof()) {
        iss.str(std::string());
        iss.clear();
        iss << line;
        getline(iss, token, ' ');
        hap.push_back(tmp);
        while (token != "-1") {
            hap.back().push_back(atoi(token.c_str()));
            getline(iss, token, ' ');
        }
        count[1].push_back(hap.back().back());
        hap.back().pop_back();
        count[0].push_back(hap.back().back());
        hap.back().pop_back();
        getline(data_file, line);
    }
    int model = floor(gsl_rng_uniform(gen) * (double)models.size());
    std::vector<double> theta, theta_prime, K, c;
    for (int i = 0; i < cfg.getRoot()["theta"].getLength(); i++) {
        theta.push_back(cfg.getRoot()["theta"][i]);
        theta_prime.push_back(cfg.getRoot()["theta_prime"][i]);
        K.push_back(cfg.getRoot()["K"][i]);
        c.push_back(cfg.getRoot()["c"][i]);
    }
    std::vector<double> theta_prior_scale = theta;
    std::vector<double> theta_prior_shape(theta.size(), 4.0);
    std::vector<double> theta_prime_prior_scale = theta_prime;
    std::vector<double> theta_prime_prior_shape(theta.size(), 4.0);
    std::vector<double> K_prior_scale = K;
    std::vector<double> K_prior_shape(K.size(), 4.0);
    std::vector<double> c_prior_scale = c;
    std::vector<double> c_prior_shape(c.size(), 4.0);
    for (unsigned int i = 0; i < theta.size(); i++) {
        theta_prior_scale[i] /= theta_prior_shape[i];
        theta_prime_prior_scale[i] /= theta_prime_prior_shape[i];
        K_prior_scale[i] /= K_prior_shape[i];
        c_prior_scale[i] /= c_prior_shape[i];
    }
    double log_w = smc(models[model], theta[model], theta_prime[model], 
        K[model], c[model], particle_numbers[model], count, hap, gen);
    std::vector<double> new_theta = theta, new_theta_prime = theta_prime, 
        new_K = K, new_c = c;
    double new_log_w = 0.0;
    int new_model = model;
    double acceptance_probability;
    int accepted = 0;
    double prop_mean = log_w;
    double prop_mean_squared = log_w * log_w;
    double sd_denom = sqrt(8.0);
    std::cout << std::setprecision(4) << std::fixed;
    for (int i = 0; i < run_length; i++) {
        new_model = floor(gsl_rng_uniform(gen) * (double)models.size());
        acceptance_probability = 0.0;
        for (unsigned int j = 0; j < theta.size(); j++) {
            new_theta[j] = theta[j] + gsl_ran_gaussian_ziggurat(gen, 
                sd / sd_denom);
            if (new_theta[j] < 0.0) {
                new_theta[j] = -new_theta[j];
            }
            acceptance_probability += log(gsl_ran_gamma_pdf(new_theta[j], 
                theta_prior_shape[j], theta_prior_scale[j])) 
                - log(gsl_ran_gamma_pdf(theta[j], theta_prior_shape[j], 
                                        theta_prior_scale[j]));
            if (models[j] > 0) {
                new_K[j] = K[j] + gsl_ran_gaussian_ziggurat(gen, sd / sd_denom);
                if (new_K[j] < 0.0) {
                    new_K[j] = -new_K[j];
                }
                acceptance_probability += log(gsl_ran_gamma_pdf(new_K[j], 
                    K_prior_shape[j], K_prior_scale[j])) 
                    - log(gsl_ran_gamma_pdf(K[j], K_prior_shape[j], 
                                            K_prior_scale[j]));
                new_c[j] = c[j] + gsl_ran_gaussian_ziggurat(gen, sd / sd_denom);
                if (new_c[j] < 0.0) {
                    new_c[j] = -new_c[j];
                }
                acceptance_probability += log(gsl_ran_gamma_pdf(new_c[j], 
                    c_prior_shape[j], c_prior_scale[j])) 
                    - log(gsl_ran_gamma_pdf(c[j], c_prior_shape[j], 
                                            c_prior_scale[j]));
                if (models[j] == 2 && mutation_flags[j] == 2) {
                    new_theta_prime[j] = theta_prime[j] 
                        + gsl_ran_gaussian_ziggurat(gen, sd / sd_denom);
                    if (new_theta_prime[j] < 0.0) {
                        new_theta_prime[j] = -new_theta_prime[j];
                    }
                    acceptance_probability 
                        += log(gsl_ran_gamma_pdf(new_theta_prime[j], 
                        theta_prime_prior_shape[j], theta_prime_prior_scale[j])) 
                        - log(gsl_ran_gamma_pdf(theta_prime[j], 
                        theta_prime_prior_shape[j], theta_prime_prior_scale[j]));
                } else if (models[j] == 2 && mutation_flags[j] == 1) {
                    new_theta_prime[j] = new_theta[j];
                } else {
                    new_theta_prime[j] = 0.0;
                }
            }
        }
        new_log_w = smc(models[new_model], theta[new_model], 
            theta_prime[new_model], K[new_model], c[new_model], 
            particle_numbers[new_model], count, hap, gen);
        prop_mean = ((double)(i + 1) * prop_mean + new_log_w) 
            / (double)(i + 2);
        prop_mean_squared = ((double)(i + 1) * prop_mean_squared + new_log_w 
            * new_log_w) / (double)(i + 2);
        acceptance_probability += new_log_w - log_w;
        if (log(gsl_rng_uniform_pos(gen)) < acceptance_probability) {
            model = new_model;
            theta = new_theta;
            theta_prime = new_theta_prime;
            K = new_K;
            c = new_c;
            log_w = new_log_w;
            accepted++;
        }
        std::cout << models[model] << "\t" << mutation_flags[model] << "\t";
        std::cout << theta[model] << "\t" << theta_prime[model] << "\t";
        std::cout << K[model] << "\t" << c[model] << "\t";
        std::cout << log_w * log10(exp(1.0)) << std::endl;
    }
    std::cout << (double)accepted / (double)run_length << std::endl;
    std::cout << prop_mean_squared - prop_mean * prop_mean << std::endl;
    return 1;
}
