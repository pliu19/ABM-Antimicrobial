#!/usr/bin/python3
# -*- coding: utf-8 -*-

import numpy as np 
import time
import json
import os 

from abm import Bacteria, Medicine
from abm import HealthCareWorker, Patient


class TestFailed(Exception):

    def __init__(self, m):
        self.message = m

    def __str__(self):
        return self.message


def draw():
    '''
    Draw a random number from 0 to 1
    '''
    return np.random.random_sample()


def bulk_initialization(args):
    '''
    Initialize all the patients at the very beginning. 
    The assignments should be dynamic based on the num patients
    '''
    patient_list = {}

    current_day = 0

    viruses = ['xa', 'xn', '0', '1', '2', '12']
    ratio_pts = [6/16, 6/16, 1/16, 1/16, 1/16, 1/16]

    patient_list[0] = Patient(args, 0, current_day, "0")
    patient_list[1] = Patient(args, 1, current_day, "1")
    patient_list[2] = Patient(args, 2, current_day, "2")
    patient_list[3] = Patient(args, 3, current_day, "12")

    for idx in range(4, args.num_patient):
        x = np.random.choice(viruses, p=ratio_pts)
        patient_list[idx] = Patient(args, idx, current_day, x)

    return patient_list


def initial_patient(args, name, current_day):
    '''
    Initialize the status of patients
    Refer to Table 0. 
    '''
    rand_num = draw()

    a = args.a 
    m = args.m 
    r1 = args.r1
    r2 = args.r2

    probs = [1-m, (1-a)*m, a*m*(1-r1-r2), a*m*r1, a*m*r2]     # Initial strain
    prCDF = [sum(probs[:idx+1]) for idx, i in enumerate(probs)]
    sList = ['xn', 'xa', '0', '1', '2']

    for pidx in range(len(prCDF)):
        if rand_num < prCDF[pidx]:
            cl = pidx
            break 
    
    patient = Patient(args, name, current_day, str(sList[cl]))
    patient.recordStatus('initialization', 0)

    return patient


def random_schedule(current_patients, num_hcw=4):
    '''
    Schedule the HCWs' routine every 24 hours  
    '''
    rst = {}

    copy_current_p = current_patients.copy()
    np.random.shuffle(copy_current_p)

    sta = 0
    end = gap = len(current_patients)//num_hcw

    for i in range(num_hcw):

        key = 'h' + str(i)
        rst[key] = current_patients[sta:end]

        sta = end
        end += gap

    return rst


def patient_colonized_a(reference_results, args, current_day, rand_num, healthW, patient, p):
    '''
    Interaction between Patient and HCW if P is colonized.
    Patient infection status = colonized and dominant strain = X_A before contact.
    Refer to Table 1.
    ToDo: Find a new prior probability.
    '''    
    hcw_strain_list = list(healthW.strain_set)

    ### hcw_strain_list may be empty
    if rand_num < p and len(hcw_strain_list) > 0:
        # Randomly choose one strain from hcw_strain_list
        chosen_strain = np.random.choice(hcw_strain_list)
        patient.dom_strain = chosen_strain
        patient.converted_infection_time(args, current_day)

        colonization_key = 'colonization_' + str(chosen_strain)
        reference_results[colonization_key] += 1
    
    else:
        # Nothing changes 
        pass 

    return healthW, patient


def patient_colonized_b(reference_results, rand_num, healthW, patient, q):
    '''
    Interaction between Patient and HCW if P is colonized.
    Patient infection status = colonized and dominant strain contains 0, 1, 2 or 12 before contact.
    Refer to Table 2.
    ToDo: Find a new prior probability.
    '''           
    ptStrain = patient.dom_strain

    if ptStrain in healthW.strain_set:
        pass 
    else:
        if rand_num < q:
            healthW.strain_set.add(ptStrain)
        else:
            pass 

    return healthW, patient


def patient_infected_a_v1(reference_results, q_rNum, r_rNum, healthW, patient, currentT, q, r):
    '''
    Interaction between Patient and HCW if P is infected.
    If status is 'I', both healthworker and patient will change their status
    Patient infection status = infected before contact and treatment time ∈ [0, 3*24] hours
    Refer to Table 3.1
    ToDo: Find a new prior probability.    
    '''
    dominant_strain = patient.dom_strain
    hcw_strain_list = list(healthW.strain_set)

    # Dominant strain is xa or xn 
    if dominant_strain == 'xa' or dominant_strain == 'xn':

        if '2' not in healthW.strain_set and '12' not in healthW.strain_set:
            pass 

        elif '2' in healthW.strain_set and '12' not in healthW.strain_set:
            if r_rNum < r:
                patient.dom_strain = '2'
                patient.super_infe = True
                patient.convt_time = currentT - patient.treat_time + 3*24
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1

        elif '12' in healthW.strain_set and '2' not in healthW.strain_set:
            if r_rNum < r:
                patient.dom_strain = '12'
                patient.super_infe = True
                patient.convt_time = currentT - patient.treat_time + 3*24
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1

        elif '2' in healthW.strain_set and '12' in healthW.strain_set: 
            if r_rNum < r:
                # equal probability to get 2 or 12 
                patient.dom_strain = np.random.choice(['2', '12'])
                patient.super_infe = True
                patient.convt_time = currentT - patient.treat_time + 3*24
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1

    else: # Dominant strain is 0, 1, 2, 12 
          # After each contact, HCW may acquire patient's dominant strain i with Probabiliy q
        if q_rNum < q:
            healthW.strain_set.add(patient.dom_strain)

        if dominant_strain == '0' or dominant_strain == '1':

            if '2' not in healthW.strain_set and '12' not in healthW.strain_set:
                pass

            elif '2' in healthW.strain_set and '12' not in healthW.strain_set:
                # Do not set conversion time when dominant_strain == 0 
                if r_rNum < r:
                    patient.dom_strain = '2'
                    if dominant_strain == '1':
                        patient.convt_time = currentT - patient.treat_time + 3*24
                        reference_results['misempiric'] += 1

                    transition_key = 'transmission_' + patient.dom_strain
                    reference_results[transition_key] += 1

            elif '12' in healthW.strain_set and '2' not in healthW.strain_set:

                if r_rNum < r:
                    patient.dom_strain = '12'
                    patient.convt_time = currentT - patient.treat_time + 3*24
                    reference_results['misempiric'] += 1
                        
                    transition_key = 'transmission_' + patient.dom_strain
                    reference_results[transition_key] += 1

            elif '2' in healthW.strain_set and '12' in healthW.strain_set:

                if r_rNum < r:
                    chosen_strain = np.random.choice(['2', '12'])
                    patient.dom_strain = chosen_strain

                    if chosen_strain == '2':
                        reference_results['tempempiric'] += 1

                        if dominant_strain == '1':
                            patient.convt_time = currentT - patient.treat_time + 3*24
                         
                    else:
                        patient.convt_time = currentT - patient.treat_time + 3*24
                        reference_results['misempiric'] += 1

                        
                    transition_key = 'transmission_' + patient.dom_strain
                    reference_results[transition_key] += 1

            else:
                raise TestFailed('Ooops')

        elif dominant_strain == '2' or dominant_strain == '12':
            pass

        else:
            raise TestFailed('Ooops')

    return healthW, patient


def patient_infected_b_v1(reference_results, q_rNum, r_rNum, healthW, patient, currentT, q, r, s):
    '''
    Interaction between Patient and HCW if P is infected.
    If status is 'I', both healthworker and patient will change their status
    Patient infection status = infected before contact and treatment time is larger than 3*24 hours
    Refer to Table 3.1
    ToDo: Find a new prior probability
    '''
    dominant_strain = patient.dom_strain
    hcw_strain_list = list(healthW.strain_set)

    # Drug use == L and dominant strain = X
    # No change for HCW, use r_rNum to determine dominant strain
    if patient.drug_use == 'L' and dominant_strain in {'xa', 'xn'}:

        if r_rNum < s and len(hcw_strain_list) > 0:
            # patient is going to change 
            chosen_strain = np.random.choice(hcw_strain_list)
            patient.dom_strain = chosen_strain

            patient.convt_time = currentT
            patient.super_infe = True
            # reference_results['superinfection'] += 1
            reference_results['cumsuperinfection'] += 1

            transition_key = 'transmission_' + patient.dom_strain
            reference_results[transition_key] += 1
            
        else:
            # no changes patient status
            pass 

    else:
        # HCW acquires patient's dominant strain i at a probability of q 
        # if not contiminated by strain i before contact
        if q_rNum < q:
            healthW.strain_set.add(patient.dom_strain)

        if patient.drug_use == 'A':
            # Under drug use = A
            if dominant_strain == '0' or dominant_strain == '2':
                # Determine which case, only three cases 
                intersection = list(healthW.strain_set.intersection({'1', '12'}))
                # HCW carries either '1' or '12'
                if len(intersection) == 1:
                    if r_rNum < r:
                        patient.dom_strain = intersection[0]
                        patient.convt_time = currentT

                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1                        

                # HCW carries both '1' and '12'
                elif len(intersection) == 2:
                    if r_rNum < r:
                        patient.dom_strain = np.random.choice(['1', '12'])
                        patient.convt_time = currentT

                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1                               

                # HCW carries neither '1' or '12'
                elif len(intersection) == 0:
                    pass 
                else:
                    raise TestFailed('Ooops')
            # Under drug use = A
            elif dominant_strain == '1' or dominant_strain == '12':
                pass
            else:
                raise TestFailed('Ooops')

        elif patient.drug_use == 'B':

            if dominant_strain == '2' or dominant_strain == '12':
                pass 

            elif dominant_strain == '0' or dominant_strain == '1':

                intersection = list(healthW.strain_set.intersection({'2', '12'}))

                if len(intersection) == 1:
                    if r_rNum < r:
                        patient.dom_strain = intersection[0]
                        patient.convt_time = currentT

                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1                          

                elif len(intersection) == 2:
                    if r_rNum < r:
                        patient.dom_strain = np.random.choice(['2', '12'])
                        patient.convt_time = currentT

                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1  

                # May conflict when HCW gets the strain from patient
                elif len(intersection) == 0:
                    pass
                else:
                    raise TestFailed('Ooops')
                    
        elif patient.drug_use == 'C':
            pass
        elif patient.drug_use == 'L' and dominant_strain not in set(['xa', 'xn']):
            pass 
        else:
            # print(patient.lab_result, dominant_strain, patient.drug_use)
            raise TestFailed('Ooops')

    return healthW, patient


def patient_infected_a_v2(reference_results, q_rNum, r_rNum, healthW, patient, currentT, q, r):
    '''
    Interaction between Patient and HCW if P is infected.
    If status is 'I', both healthworker and patient will change their status
    Patient infection status = infected before contact and treatment time ∈ [0, 3*24] hours
    Refer to Table 3.1
    ToDo: Find a new prior probability.    
    '''
    dominant_strain = patient.dom_strain
    hcw_strain_list = list(healthW.strain_set)

    # Dominant strain is xa or xn 
    if dominant_strain == 'xa' or dominant_strain == 'xn':

        if '2' not in healthW.strain_set and '12' not in healthW.strain_set:
            pass 

        elif '2' in healthW.strain_set and '12' not in healthW.strain_set:
            if r_rNum < r:
                patient.dom_strain = '2'
                patient.super_infe = True
                patient.convt_time = currentT - patient.treat_time + 3*24
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1  

        elif '12' in healthW.strain_set and '2' not in healthW.strain_set:
            if r_rNum < r:
                patient.dom_strain = '12'
                patient.super_infe = True
                patient.convt_time = currentT - patient.treat_time + 3*24
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1  

        elif '2' in healthW.strain_set and '12' in healthW.strain_set: 
            if r_rNum < r:
                # equal probability to get 2 or 12 
                patient.dom_strain = np.random.choice(['2', '12'])
                patient.super_infe = True
                patient.convt_time = currentT - patient.treat_time + 3*24
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1  

    else: # Dominant strain is 0, 1, 2, 12 
          # After each contact, HCW may acquire patient's dominant strain i with Probabiliy q
        if q_rNum < q:
            healthW.strain_set.add(dominant_strain)

        if dominant_strain == '0' or dominant_strain == '1':

            if '2' not in healthW.strain_set and '12' not in healthW.strain_set:
                pass

            elif '2' in healthW.strain_set and '12' not in healthW.strain_set:
                # Do not set conversion time when dominant_strain == 0 
                # This is not true for v2 
                if r_rNum < r:
                    patient.dom_strain = '2'
                    patient.convt_time = currentT - patient.treat_time + 3*24
                    reference_results['misempiric'] += 1
                    transition_key = 'transmission_' + patient.dom_strain
                    reference_results[transition_key] += 1  

            elif '12' in healthW.strain_set and '2' not in healthW.strain_set:

                if r_rNum < r:
                    patient.dom_strain = '12'
                    patient.convt_time = currentT - patient.treat_time + 3*24
                    reference_results['misempiric'] += 1
                    transition_key = 'transmission_' + patient.dom_strain
                    reference_results[transition_key] += 1  

            elif '2' in healthW.strain_set and '12' in healthW.strain_set:

                if r_rNum < r:
                    chosen_strain = np.random.choice(['2', '12'])
                    patient.dom_strain = chosen_strain
                    patient.convt_time = currentT - patient.treat_time + 3*24
                    transition_key = 'transmission_' + patient.dom_strain
                    reference_results[transition_key] += 1  

            else:
                raise TestFailed('Ooops')

        elif dominant_strain == '2' or dominant_strain == '12':
            pass

        else:
            raise TestFailed('Ooops')

    return healthW, patient


def patient_infected_b_v2(reference_results, q_rNum, r_rNum, healthW, patient, currentT, q, r, s):
    '''
    Interaction between Patient and HCW if P is infected.
    If status is 'I', both healthworker and patient will change their status
    Patient infection status = infected before contact and treatment time is larger than 3*24 hours
    Refer to Table 3.1
    ToDo: Find a new prior probability
    '''
    dominant_strain = patient.dom_strain
    hcw_strain_list = list(healthW.strain_set)

    # Drug use == L and dominant strain = X
    # No change for HCW, use r_rNum to determine dominant strain
    if patient.drug_use == 'B' and dominant_strain in {'xa', 'xn'}:

        if '2' not in healthW.strain_set and '12' not in healthW.strain_set:
            pass 

        elif '2' in healthW.strain_set and '12' not in healthW.strain_set:
            if r_rNum < r:
                patient.dom_strain = '2'
                patient.super_infe = True
                patient.convt_time = currentT
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                # reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1  

        elif '12' in healthW.strain_set and '2' not in healthW.strain_set:
            if r_rNum < r:
                patient.dom_strain = '12'
                patient.super_infe = True
                patient.convt_time = currentT
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                # reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1  

        elif '2' in healthW.strain_set and '12' in healthW.strain_set: 
            if r_rNum < r:
                # equal probability to get 2 or 12 
                patient.dom_strain = np.random.choice(['2', '12'])
                patient.super_infe = True
                patient.convt_time = currentT
                # reference_results['superinfection'] += 1
                reference_results['cumsuperinfection'] += 1
                # reference_results['misempiric'] += 1
                transition_key = 'transmission_' + patient.dom_strain
                reference_results[transition_key] += 1  

        # if r_rNum < s and len(hcw_strain_list) > 0:
        #     # patient is going to change 
        #     chosen_strain = np.random.choice(hcw_strain_list)
        #     patient.dom_strain = chosen_strain

        #     patient.convt_time = currentT
        #     patient.super_infe = True
        #     reference_results['superinfection'] += 1
        #     reference_results['cumsuperinfection'] += 1
        #     transition_key = 'transmission_' + patient.dom_strain
        #     reference_results[transition_key] += 1  

        # else:
        #     # no changes patient status
        #     pass 

    else:
        # HCW acquires patient's dominant strain i at a probability of q 
        # if not contiminated by strain i before contact
        if q_rNum < q:
            healthW.strain_set.add(patient.dom_strain)

        if patient.drug_use == 'A':
            # Under drug use = A
            if dominant_strain == '0' or dominant_strain == '2':
                # Determine which case, only three cases 
                intersection = list(healthW.strain_set.intersection({'1', '12'}))
                # HCW carries either '1' or '12'
                if len(intersection) == 1:
                    if r_rNum < r:
                        patient.dom_strain = intersection[0]
                        patient.convt_time = currentT
                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1  

                # HCW carries both '1' and '12'
                elif len(intersection) == 2:
                    if r_rNum < r:
                        patient.dom_strain = np.random.choice(['1', '12'])
                        patient.convt_time = currentT
                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1  

                # HCW carries neither '1' or '12'
                elif len(intersection) == 0:
                    pass 
                else:
                    raise TestFailed('Ooops')
            # Under drug use = A
            elif dominant_strain == '1' or dominant_strain == '12':
                pass
            else:
                raise TestFailed('Ooops')

        elif patient.drug_use == 'B':

            if dominant_strain == '2' or dominant_strain == '12':
                pass 

            elif dominant_strain == '0' or dominant_strain == '1':

                intersection = list(healthW.strain_set.intersection({'2', '12'}))

                if len(intersection) == 1:
                    if r_rNum < r:
                        patient.dom_strain = intersection[0]
                        patient.convt_time = currentT
                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1  

                elif len(intersection) == 2:
                    if r_rNum < r:
                        patient.dom_strain = np.random.choice(['2', '12'])
                        patient.convt_time = currentT
                        transition_key = 'transmission_' + patient.dom_strain
                        reference_results[transition_key] += 1  

                # May conflict when HCW gets the strain from patient
                elif len(intersection) == 0:
                    pass
                else:
                    raise TestFailed('Ooops')
                    
        elif patient.drug_use == 'C':
            pass

        else:
            raise TestFailed('Ooops')

    return healthW, patient



def infection_development(patients_record, current_patients, current_day, reference_results):
    '''
    Patient may convert from 'C' to 'I' that was determined when initializing the patient. 
    '''
    for key in current_patients:
        patient = patients_record[key]
        if patient.status == 'C':

            if patient.infct_flag and patient.infct_time == current_day:

                patients_record[key].status = 'I'
                patients_record[key].treat_time = 0
                patients_record[key].drug_use = 'B'
                if patient.dom_strain in {'0', '1', '2', '12'}:
                    patients_record[key].lab_result = patients_record[key].dom_strain
                    name2 = 'cuminfection_' + patient.dom_strain
                    reference_results[name2] += 1
                    key_name = 'labresult_' + str(patient.dom_strain)
                    reference_results[key_name] += 1
                else:
                    patients_record[key].lab_result = 'X'
                    reference_results['cuminfection_X'] += 1


    return patients_record


def hcw_cleanUp(healthW, currentTime, eta=0.5):
    '''
    Health worker's strains may be cleaned up
    ToDo: Set eta to clean up all the strains from HCW 
    '''
    rand_num = draw()
    if currentTime%8 == 0:

        healthW.strain_set.clear()
        return healthW
    
    if rand_num < eta:
        healthW.strain_set.clear()

    return healthW


def increase_treatment_icu_time(patients_record, current_patients):
    '''
    Increase the ICU time and treatment time if exists. 
    '''
    for key in current_patients:
        patient = patients_record[key]
        patient.time_inICU += 24
        
        # Wrong, isinstance
        if isinstance(patient.treat_time, int):
            patients_record[key].treat_time += 24

    return patients_record    


def drug_change_v1(reference_results, patients_record, current_patients, currentTime):
    '''
    Drug in use may change under some conditions. 
    '''
    for key in current_patients:
        patient = patients_record[key]
        if patient.treat_time == 3*24:
            
            if patient.lab_result == '0' or patient.lab_result == '2':
                patients_record[key].drug_use = 'A'
                reference_results['def_drug_use_A'] += 1
            
            elif patient.lab_result == '1':
                patients_record[key].drug_use = 'B'
                reference_results['def_drug_use_B'] += 1

            elif patient.lab_result == '12':
                patients_record[key].drug_use = 'C'
                reference_results['def_drug_use_C'] += 1

            elif patient.lab_result == 'X':
                patients_record[key].drug_use = 'L'
                reference_results['def_drug_use_L'] += 1

            else:
                raise TestFailed('Ooops')

        if patient.convt_time and (currentTime - patient.convt_time) == 3*24:

            if patient.drug_use == 'A' or patient.drug_use == 'B':
                patients_record[key].drug_use = 'C'
                reference_results['corr_drug_use_C'] += 1
            elif patient.drug_use == 'L':
                patients_record[key].drug_use = 'B'
                patients_record[key].treat_time = 0
                patients_record[key].lab_result = patients_record[key].dom_strain
                reference_results['corr_drug_use_B'] += 1

    return patients_record


def drug_change_v2(reference_results, patients_record, current_patients, currentTime):
    '''
    Drug in use may change under some conditions. 
    '''
    for key in current_patients:
        patient = patients_record[key]
        if patient.treat_time == 3*24:
            
            if  patient.lab_result == '2':
                patients_record[key].drug_use = 'A'
                reference_results['def_drug_use_A'] += 1
            
            elif patient.lab_result == '1' or patient.lab_result == '0' or patient.lab_result == 'X':
                patients_record[key].drug_use = 'B'
                reference_results['def_drug_use_B'] += 1

            elif patient.lab_result == '12':
                patients_record[key].drug_use = 'C'
                reference_results['def_drug_use_C'] += 1

            else:
                raise TestFailed('Ooops')

        if patient.convt_time and (currentTime - patient.convt_time) == 3*24:

            patients_record[key].drug_use = 'C'
            reference_results['corr_drug_use_C'] += 1

    return patients_record


def check_treatment_completion(patient, currentT):
    '''
    Check whether the treatment of patient is complete 
    '''
    if patient.convt_time is None and \
        patient.lab_result in ['0', '1'] and \
        patient.treat_time == 7*24:

        return True

    if patient.convt_time is None and \
        patient.lab_result in ['2', '12'] and \
        patient.treat_time == 10*24:

        return True

    if patient.convt_time is not None and \
        currentT - patient.convt_time == 10*24: 

        return True

    return False


def treatment_completion(reference_results, patients_record, current_patients, currentT):

    for key in current_patients:
        patient = patients_record[key]

        if check_treatment_completion(patient, currentT):
            # This part is the decrement part of superinfection
            # super_infection = patients_record[key].super_infe
            # if super_infection:
            #     reference_results['superinfection'] -= 1

            patients_record[key].status = 'C' 
            patients_record[key].treat_time = None
            patients_record[key].convt_time = None
            patients_record[key].lab_result = None
            patients_record[key].drug_use = None

            # patients_record[key].recordStatus('Treatment Completion', currentT)
        else:
            pass 

    return patients_record


def discharge_admission(patients_record, currentT, args, discharge_probas, current_patients, patient_last_idx, reference_results):
    '''
    Discharged happens when
    Has infection status=colonized
    Has time in ICU ≥ 2(days)
    '''
    current_day = currentT // 24

    kappa_mu = args.kappa_mu
    copy_current_patients = current_patients.copy()

    for key in copy_current_patients:

        patient = patients_record[key]

        if patient.time_inICU == 40*24:
            # Admission of new patient
            patient_last_idx += 1
            new_patient = initial_patient(args, patient_last_idx, current_day)
            current_patients.remove(key)
            current_patients.append(patient_last_idx)
            patients_record[patient_last_idx] = new_patient
            reference_results['admission'] += 1
            reference_results['discharge'] += 1
            # super infection discount 
            # super_infection = patients_record[key].super_infe
            # if super_infection:
            #     reference_results['superinfection'] -= 1

            continue

        if patient.super_infe:
            continue

        if patient.status is 'C' and patient.time_inICU >= 2*24:

            day = patient.time_inICU // 24
            rNum = draw()
            if rNum < discharge_probas[day]:
                # Admission of new patient
                patient_last_idx += 1
                new_patient = initial_patient(args, patient_last_idx, current_day)
                current_patients.remove(key)
                current_patients.append(patient_last_idx)
                patients_record[patient_last_idx] = new_patient
                reference_results['admission'] += 1
                reference_results['discharge'] += 1
                # super infection discount 
                # super_infection = patients_record[key].super_infe
                # if super_infection:
                #     reference_results['superinfection'] -= 1

        elif patient.status is 'I' and patient.time_inICU >= 2*24:

            day = patient.time_inICU // 24
            rNum = draw()
            if rNum < kappa_mu * discharge_probas[day]:
                # Admission of new patient
                patient_last_idx += 1
                new_patient = initial_patient(args, patient_last_idx, current_day)
                current_patients.remove(key)
                current_patients.append(patient_last_idx)
                patients_record[patient_last_idx] = new_patient
                reference_results['admission'] += 1
                reference_results['discharge'] += 1
                # super infection discount 
                # super_infection = patients_record[key].super_infe
                # if super_infection:
                #     reference_results['superinfection'] -= 1

    return patients_record, current_patients, patient_last_idx


def death_event(patients_record, current_patients, patient_last_idx, currentT, args, death_probs, reference_results):
    '''
    Death occurs at a probability of v1, v2, v3 for three situations. 
    Obviously, v1 < v2 < v3
    '''
    current_day = currentT // 24
    kappa_nu = args.kappa_nu

    for key in current_patients:

        patient = patients_record[key]
        day = patient.time_inICU // 24

        if patient.convt_time is None or currentT - patient.convt_time > 3*24:

            if draw() < death_probs[day]:
                # Admission of new patient
                patient_last_idx += 1
                new_patient = initial_patient(args, patient_last_idx, current_day)
                current_patients.remove(key)
                current_patients.append(patient_last_idx)
                patients_record[patient_last_idx] = new_patient
                reference_results['admission'] += 1
                reference_results['death'] += 1
                # super infection discount 
                # super_infection = patients_record[key].super_infe
                # if super_infection:
                #     reference_results['superinfection'] -= 1
                
        elif patient.convt_time is not None and currentT - patient.convt_time <= 3*24:
            
            if draw() < death_probs[day] * kappa_nu:
                # Admission of new patient
                patient_last_idx += 1
                new_patient = initial_patient(args, patient_last_idx, current_day)
                current_patients.remove(key)
                current_patients.append(patient_last_idx)
                patients_record[patient_last_idx] = new_patient
                reference_results['admission'] += 1
                reference_results['death'] += 1
                # super infection discount 
                # super_infection = patients_record[key].super_infe
                # if super_infection:
                #     reference_results['superinfection'] -= 1

        else:
            raise ('Error')

    return patients_record, current_patients, patient_last_idx


def intrinsic_mutation_v1(patients_record, current_patients, args, currentT, reference_results):
    '''
    Intrinsic mutation process
    '''
    epsilon = args.epsilon
    copy_current_patients = current_patients.copy()

    for key in copy_current_patients:

        rand_num = draw() 
        patient = patients_record[key]

        if patient.status == 'I':

            if patient.treat_time <= 3*24:
                if rand_num < epsilon:
                    if patient.dom_strain == '0':
                        patients_record[key].dom_strain = '2'
                        reference_results['mutation_2'] += 1
                    elif patient.dom_strain == '1':
                        patients_record[key].dom_strain = '12'
                        patients_record[key].convt_time = currentT - patient.treat_time + 3*24
                        reference_results['mutation_12'] += 1
                    # elif patient.dom_strain in {'xa', 'xn'}:
                    #     patients_record[key].dom_strain = '2'
                    #     patients_record[key].convt_time = currentT - patient.treat_time + 3*24
                    #     patients_record[key].super_infe = True
                    #     reference_results['superinfection'] += 1
                    #     reference_results['mutation_2'] += 1
            else:
                if rand_num < epsilon:
                    if patient.dom_strain == '0' and patient.drug_use == 'A':
                        patients_record[key].dom_strain = '1'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_1'] += 1
                    elif patient.dom_strain == '2' and patient.drug_use == 'A':
                        patients_record[key].dom_strain = '12'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_12'] += 1
                    elif patient.dom_strain == '0' and patient.drug_use == 'B':
                        patients_record[key].dom_strain = '2'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_2'] += 1
                    elif patient.dom_strain == '1' and patient.drug_use == 'B':
                        patients_record[key].dom_strain = '12'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_12'] += 1

    return patients_record 


def intrinsic_mutation_v2(patients_record, current_patients, args, currentT, reference_results):
    '''
    Intrinsic mutation process
    '''
    epsilon = args.epsilon
    copy_current_patients = current_patients.copy()

    for key in copy_current_patients:

        rand_num = draw() 
        patient = patients_record[key]

        if patient.status == 'I':

            if patient.treat_time <= 3*24:
                if rand_num < epsilon:
                    if patient.dom_strain == '0':
                        patients_record[key].dom_strain = '2'
                        patients_record[key].convt_time = currentT - patient.treat_time + 3*24  # The difference from v1
                        reference_results['mutation_2'] += 1
                    elif patient.dom_strain == '1':
                        patients_record[key].dom_strain = '12'
                        patients_record[key].convt_time = currentT - patient.treat_time + 3*24
                        reference_results['mutation_12'] += 1
                    # elif patient.dom_strain in {'xa', 'xn'}:
                    #     patients_record[key].dom_strain = '2'
                    #     patients_record[key].convt_time = currentT - patient.treat_time + 3*24
                    #     patients_record[key].super_infe = True
                    #     reference_results['superinfection'] += 1
                    #     reference_results['mutation_2'] += 1
            else:
                if rand_num < epsilon:
                    if patient.dom_strain == '0' and patient.drug_use == 'A':
                        patients_record[key].dom_strain = '1'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_1'] += 1
                    elif patient.dom_strain == '2' and patient.drug_use == 'A':
                        patients_record[key].dom_strain = '12'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_12'] += 1
                    elif patient.dom_strain == '0' and patient.drug_use == 'B':
                        patients_record[key].dom_strain = '2'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_2'] += 1
                    elif patient.dom_strain == '1' and patient.drug_use == 'B':
                        patients_record[key].dom_strain = '12'
                        patients_record[key].convt_time = currentT
                        reference_results['mutation_12'] += 1
                    # elif patient.dom_strain in {'xa', 'xn'} and patient.drug_use == 'B':
                    #     patients_record[key].dom_strain = '2'
                    #     patients_record[key].convt_time = currentT
                    #     patients_record[key].super_infe = True

    return patients_record 

############### Need to take another look and rewrite under this line


