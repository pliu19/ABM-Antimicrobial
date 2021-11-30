#!/usr/bin/python3
# -*- coding: utf-8 -*-
import os
import copy
import argparse

import math
import time
import pickle
from itertools import product

import scipy.integrate as integrate
import scipy.stats as stats

from abm import Bacteria, Medicine
from abm import HealthCareWorker, Patient
from scheme_func import draw, initial_patient, random_schedule, bulk_initialization
from scheme_func import patient_colonized_a, patient_colonized_b
from scheme_func import patient_infected_a_v1, patient_infected_b_v1
from scheme_func import patient_infected_a_v2, patient_infected_b_v2
from scheme_func import infection_development, drug_change_v1, hcw_cleanUp
from scheme_func import increase_treatment_icu_time, treatment_completion
from scheme_func import discharge_admission, death_event, intrinsic_mutation_v1
from scheme_func import drug_change_v2, intrinsic_mutation_v2


def truncnorm_func(mean, std=0.01):

    min_epsilon = 10e-5

    bounded_length = mean * 0.1

    lower, upper = mean-bounded_length, mean+bounded_length

    X = stats.truncnorm((lower - mean) / std, (upper - mean) / std, loc=mean, scale=std)

    return max(min_epsilon, X.rvs())


def function(x, lambdat):
    
    return lambdat


def read_csv(file_path):
    '''
    Read the csv files to get probabilities
    '''
    record = {}
    result = {}

    with open(file_path, 'r') as f:
        lis = [line.split() for line in f]        
        for i, x in enumerate(lis):              
            day_idx, prob = x[0].split(',')
            day_idx = int(day_idx)
            prob = float(prob)
            record[day_idx] = prob
    
    for key, value in record.items():

        I = integrate.quad(function, 0, key, args=(value))[0]
        proba = value * math.exp(-I)
        result[key] = proba

    return result


def ph_interaction(patient, healthW, args, current_day, current_hour, reference_results, func_map):
    """
    patient healthcare interaction function
    """
    strain = patient.dom_strain
    status = patient.status
    trtime = patient.treat_time

    p = args.p
    q = args.q
    r = args.r
    s = args.s 

    if status == "C" and strain == "xa":

        rand_num = draw()
        healthW, patient = patient_colonized_a(reference_results, args, current_day, rand_num, healthW, patient, p=p)

    elif status == "C" and strain in {"0", "1", "2", "12"}:

        rand_num = draw()
        healthW, patient = patient_colonized_b(reference_results, rand_num, healthW, patient, q=q)

    elif status == "I" and trtime and trtime <= 3*24:

        q_rNum = draw()
        r_rNum = draw()
        healthW, patient = func_map['patient_infected_a'](reference_results, q_rNum, r_rNum, healthW, patient, current_hour, q=q, r=r)

    elif status == "I" and trtime and trtime > 3*24:

        q_rNum = draw()
        r_rNum = draw()
        healthW, patient = func_map['patient_infected_b'](reference_results, q_rNum, r_rNum, healthW, patient, current_hour, q=q, r=r, s=s)

    else:
        pass

    return healthW, patient


def recordEverything(file, version, file_path='./log/'):

    timestr = time.strftime('%Y%m%d-%H%M%S')
    fileName = file_path + '/' + version  + '_' + timestr + '.pkl'

    pickle.dump(file, open(fileName, "wb" ) )

def main(args, func_map):

    reference_results = {}
    reference_results['admission'] = 0         # Done
    reference_results['death'] = 0             # Done
    reference_results['discharge'] = 0         # Done
    reference_results['cuminfection_0'] = 0    # Done  
    reference_results['cuminfection_1'] = 0    # Done
    reference_results['cuminfection_2'] = 0    # Done
    reference_results['cuminfection_12'] = 0   # Done
    reference_results['cuminfection_X'] = 0    # Done 
    reference_results['labresult_0'] = 0       # Done. only happens in infection_development
    reference_results['labresult_1'] = 0       # Done. only happens in infection_development
    reference_results['labresult_2'] = 0       # Done. only happens in infection_development
    reference_results['labresult_12'] = 0      # Done. only happens in infection_development
    # reference_results['superinfection'] = 0    # Done. 
    reference_results['cumsuperinfection'] = 0 # Done 
    reference_results['colonization_0'] = 0    # Done
    reference_results['colonization_1'] = 0    # Done
    reference_results['colonization_2'] = 0    # Done
    reference_results['colonization_12'] = 0   # Done
    reference_results['misempiric'] = 0        # Done
    reference_results['tempempiric'] = 0       # Done. Only one condition. 

    reference_results['def_drug_use_A'] = 0        # Done
    reference_results['def_drug_use_B'] = 0        # Done
    reference_results['def_drug_use_C'] = 0        # Done
    reference_results['def_drug_use_L'] = 0        # Done

    reference_results['corr_drug_use_A'] = 0        # Done
    reference_results['corr_drug_use_B'] = 0        # Done
    reference_results['corr_drug_use_C'] = 0        # Done
    reference_results['corr_drug_use_L'] = 0        # Done

    reference_results['mutation_1']  = 0
    reference_results['mutation_2']  = 0
    reference_results['mutation_12']  = 0
    reference_results['transmission_0'] = 0
    reference_results['transmission_1'] = 0
    reference_results['transmission_2'] = 0
    reference_results['transmission_12'] = 0

    final_results = {}
    final_results['infection_0'] = [0] * args.days
    final_results['infection_1'] = [0] * args.days
    final_results['infection_2'] = [0] * args.days
    final_results['infection_12'] = [0] * args.days
    final_results['infection_X'] = [0] * args.days
    final_results['superinfection'] = [0] * args.days

    for key, value in reference_results.items():
        final_results[key] = []

    current_hour = 0 
    current_day = 0 

    num_day = args.days
    num_pat = args.num_patient
    num_hcw = args.num_hcw

    hcw_list = {}

    death_probs = read_csv(args.death_prob_file)
    discharge_probs = read_csv(args.discharge_prob_file)

    assert args.num_patient % args.num_hcw == 0

    # Initialize patients 
    patient_list = bulk_initialization(args)

    patient_idx = args.num_patient - 1
    current_p_names = list(range(num_pat))

    # Initialize hcw 
    for i in range(num_hcw):
        name = 'h' + str(i)
        hcw_list[name] = HealthCareWorker()

    # Simulation event
    for _ in range(num_day):
        
        # Every 8 hour is a shift
        for shift in range(3):
            # The keys of schedule are h0, h1, h2, h3...
            schedule = random_schedule(current_p_names, num_hcw)
            visit_idx = 0 
            
            # For example, 16 patients and 4 hcws. 
            for _ in range(args.num_patient // args.num_hcw):  
                # Should operate at the same time
                for hcw_idx in range(args.num_hcw):

                    hcw_id = f"h{hcw_idx}"
                    h_pt_idx = schedule[hcw_id][visit_idx]
                    h_patient = patient_list[h_pt_idx]

                    hcw_list[hcw_id], patient_list[h_pt_idx] = ph_interaction(h_patient, hcw_list[hcw_id], args, current_day, current_hour, reference_results, func_map)

                # Update the results and record 
                visit_idx += 1 
                # Increase the time interval by time_interval (e.g., 2)
                current_hour += args.time_interval
        
            # 6. HCW clean up the virus
            for key, value in hcw_list.items():
                hcw_list[key] = hcw_cleanUp(value, current_hour, args.eta)

        # Other events which happen once per day

        # 1. infection development
        patient_list = infection_development(patient_list, current_p_names, current_day, reference_results)

        # 2. intrinsic mutation
        patient_list = func_map['intrinsic_mutation'](patient_list, current_p_names, args, current_hour, reference_results)

        # 3. Increase the icu time and treatment time, since it is the end of day. Should be updated first
        patient_list = increase_treatment_icu_time(patient_list, current_p_names)

        # 4. drug change 
        patient_list = func_map['drug_change'](reference_results, patient_list, current_p_names, current_hour)

        # 5. completion of treatment 
        patient_list = treatment_completion(reference_results, patient_list, current_p_names, current_hour)

        # 6. discharge and admission
        patient_list, current_p_names, patient_idx = discharge_admission(patient_list, current_hour, args, discharge_probs, current_p_names, patient_idx, reference_results)

        # 7. deaths
        patient_list, current_p_names, patient_idx = death_event(patient_list, current_p_names, patient_idx, current_hour, args, death_probs, reference_results)         

        for key, value in reference_results.items():
            final_results[key].append(value)

        for key in current_p_names:
            patient = patient_list[key]
            # Check for patient's infenction
            if patient.status == 'I':
                strain = 'X' if patient.dom_strain in {'xa', 'xn'} else patient.dom_strain 
                infection_name = 'infection_' + strain
                final_results[infection_name][current_day] += 1

            # Check for patient's super infection
            if patient.super_infe:
            	final_results['superinfection'][current_day] += 1


        current_day += 1

    return final_results

     
if __name__ == "__main__":


    parser = argparse.ArgumentParser(description='Parameters of simulation script')

    parser.add_argument('--version', default='v1',
                        help='the version of simulation')

    parser.add_argument('--days', type=int, default=1460,
                        help='the days of simulation')

    parser.add_argument('--time_interval', default=2,
                        help='the interval time that a HCW visit for a patient')

    parser.add_argument('--num_patient', type=int, default=16,
                        help='the number of patients, this is fixed for now. Do not change!')

    parser.add_argument('--num_hcw', type=int, default=4,
                        help='the number of patients, this is fixed for now. Do not change!')

    parser.add_argument('--num_strain', default=6,
                        help='the number of strains (0, 1, 2, 12, xa, xn)')

    parser.add_argument('--num_drug', default=4,
                        help='the number of drugs (A, B, C, L)')

    parser.add_argument('--std', default=0.0,
                        help='the standard deviation of all the parameters. If 0, then the parameter is deterministic')

    parser.add_argument('--a', default=0.1,
                        help='the parameter to initialize a patient, the rest of probas should be (1-ma)')

    parser.add_argument('--m', default=0.6,
                        help='the parameter to initialize a patient')

    parser.add_argument('--r1', default=0.35,
                        help='the parameter to initialize a patient')    

    parser.add_argument('--r2', default=0.23,
                        help='the parameter to initialize a patient')

    parser.add_argument('--p', default=0.1,
                        help='the probability of colonization for uninfected patients')

    parser.add_argument('--q', default=0.05,
                        help='the probability of contamination for HCW in each contact')

    parser.add_argument('--r', default=0.1,
                        help='the probability of colonization for infected patients')

    parser.add_argument('--s', default=0.015,
                        help='the probability of super-infection')

    parser.add_argument('--sigmax', default=0.16,
                        help='the probability for infection development of strain X')

    parser.add_argument('--sigmac', default=0.45,
                        help='the probability for infection development of ARB strains')

    parser.add_argument('--epsilon', default=0.03,
                        help='the probability of mutation')

    parser.add_argument('--eta', default=0.5,
                        help='the probability of HCW to clean up the strains after every contact')

    parser.add_argument('--discharge_prob_file', default='../dischargeprob.csv',
                        help='the file that contains the Probability of discharge for colonized patients each day')

    parser.add_argument('--kappa_mu', default=0.74,
                        help='hazard ratio of discharge for infected patients')

    parser.add_argument('--death_prob_file', default='../deathprob.csv',
                        help='the probability of death for colonized patients')

    parser.add_argument('--kappa_nu', default=1.04,
                        help='hazard ratio of death for patients with inadequate antimicrobials')

    parser.add_argument('--num_runs', default=1000, type=int, 
                        help='number of runs for the simulation')

    args = parser.parse_args()

    func_map = {}

    if args.version == 'v1':
        func_map['drug_change'] = drug_change_v1
        func_map['patient_infected_a'] = patient_infected_a_v1
        func_map['patient_infected_b'] = patient_infected_b_v1
        func_map['intrinsic_mutation'] = intrinsic_mutation_v1
    elif args.version == 'v2':
        func_map['drug_change'] = drug_change_v2
        func_map['patient_infected_a'] = patient_infected_a_v2
        func_map['patient_infected_b'] = patient_infected_b_v2
        func_map['intrinsic_mutation'] = intrinsic_mutation_v2    
    else:
        raise ("No such option")

    original_args = copy.deepcopy(args)

    p_q_list = [(0.1, 0.1), (0.05, 0.3), (0.3, 0.05), (0.05, 0.05)]

    p_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
    q_list = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]

    for p, q in list(product(p_list, q_list)):

    # for p, q in p_q_list:

        args.q = q
        args.p = p 
        args.r = p 

        whole_res = {}
        file_path = './log/' + args.version + '_q=' + str(args.q)+ '_p&r=' + str(args.p) + '/'

        if not os.path.exists(file_path):
            os.makedirs(file_path)

        for i in range(args.num_runs):
            
            if args.std > 0.:
                args.a = truncnorm_func(original_args.a, args.std)
                args.m = truncnorm_func(original_args.m, args.std)
                args.r1 = truncnorm_func(original_args.r1, args.std)
                args.r2 = truncnorm_func(original_args.r2, args.std)
                args.p = truncnorm_func(original_args.p, args.std)
                args.q = truncnorm_func(original_args.q, args.std)
                args.r = truncnorm_func(original_args.r, args.std)
                args.s = truncnorm_func(original_args.s, args.std)
                args.epsilon = truncnorm_func(original_args.epsilon, args.std)
                args.sigmax = truncnorm_func(original_args.sigmax, args.std)
                args.sigmac = truncnorm_func(original_args.sigmac, args.std)
                args.eta = truncnorm_func(original_args.eta, args.std)
                args.kappa_mu = truncnorm_func(original_args.kappa_mu, args.std)
                args.kappa_nu = truncnorm_func(original_args.kappa_nu, args.std)

            res = main(args, func_map)
            whole_res[i] = res

            if i % 10 == 0:
                print(p, q, i)

        recordEverything(whole_res, args.version, file_path)
