#!/usr/bin/python3
# -*- coding: utf-8 -*-
import random
import numpy as np

def draw():
    '''
    Draw a random number from 0 to 1
    '''
    return np.random.random_sample()


class Bacteria:
    '''
    This is the class for Bacteria.
    '''
    def __init__(self, name):
        self.name = name 


class Medicine:
    '''
    This is the class for Medicine.
    '''
    def __init__(self, name):
        self.name = name 


class HealthCareWorker:
    '''
    This is the class for Health Care Worker (HCW).

    Attributes:
        strain (Set(String)): The set of strans a HCW carries
        record:  The record of HCW from beginning to the end at each timestamp
    '''

    def __init__(self, strain_set=None):
        '''
        Initialize a HCW given the strain set or not. 

        Parameters:
            strain_set: set object
            record: All the records during the simulation
        Returns:
            None
        '''
        if strain_set:
            self.strain_set = strain_set
        else:
            self.strain_set = set()
        self.record = {}

    def recordStatus(self, reason, ts):
        '''
        Record the status at each time. 
        
        Parameters:
            Reason: a text meassage to record the record of status changing. 
            ts: the time stamp 
        Returns:
            None
        '''        
        if len(self.record) == 0:
            idx = 0
        else:
            idx = max(self.record.keys()) + 1

        self.record[idx] = {}
        self.record[idx]['strain_set'] = sorted(list(self.strainSet))
        self.record[idx]['ts'] = ts
        self.record[idx]['reason'] = reason


class Patient:
    """
    This is the class for Patient.

    Attributes:
        status (str): could be "C" (colonized) or "I" (infected), initialed by "C"
        dom_strain (str): 0, 1, 2, 12, X, (X includes X_A, X_N)
        drug_use (str): 
        treat_time (int): increases w.r.t. time once set to be 0, initialed by None
        convt_time (int): does not increase w.r.t. time, initialed by None
        lab_result (str): 0, 1, 2, 12, X or Null, initialed by None
        time_inICU (int): increases w.r.t. time, initialed by 0 
        super_infe (bool): True or False, initial by False
        record: All the records during the simulation
    """
    def __init__(self, args, name, current_day, strain=None):
        """
        Parameters:
            strain (String): to initial the value of self.dom_strain
        Returns:
            None: None
        """
        self.name = name
        self.status = "C"
        self.dom_strain = strain
        self.drug_use = None
        self.infct_time = None    # in days
        self.infct_flag = False   # Not going to be infected by default
        self.treat_time = None
        self.convt_time = None
        self.lab_result = None 
        self.time_inICU = 0
        self.super_infe = False 

        self.record = {}

        self.infection_time(args, current_day)


    def infection_time(self, args, current_day, time_interval=5):
        '''
        Infection development
        '''

        if self.dom_strain == 'xa' or self.dom_strain == 'xn':
            sigma = args.sigmax
        else:
            sigma = args.sigmac

        rand_num = draw()

        if rand_num < sigma:
            # Going to be infected
            self.infct_flag = True
            self.infct_time = current_day + random.choice(list(range(0, time_interval+1))) 
        else:
            # Not going to be infected 
            pass 

    
    def converted_infection_time(self, args, current_day, time_interval=5):
        '''
        If patient convert from X to others, re-determine the infection time
        '''
        sigma = args.sigmac

        rand_num = draw()

        if rand_num < sigma:
            # Going to be infected
            self.infct_flag = True
            self.infct_time = random.choice(list(range(0, time_interval+1))) + current_day 
        else:
            # Not going to be infected 
            pass 


    def recordStatus(self, reason, currentT):
        '''
        Record the status at each time 
        '''
        if len(self.record) == 0:
            idx = 0
        else:
            idx = max(self.record.keys()) + 1

        self.record[idx] = {}
        self.record[idx]["attributes"] = {}
        self.record[idx]['reason'] = reason
        self.record[idx]['currentTS'] = currentT

        self.record[idx]["attributes"]["name"] = self.name
        self.record[idx]["attributes"]["status"] = self.status
        self.record[idx]["attributes"]["drug_use"] = self.drug_use
        self.record[idx]["attributes"]["dom_strain"] = self.dom_strain
        self.record[idx]["attributes"]["lab_result"] = self.lab_result
        self.record[idx]["attributes"]["time_inICU"] = self.time_inICU
        self.record[idx]["attributes"]["infct_time"] = self.infct_time
        self.record[idx]["attributes"]["treat_time"] = self.treat_time
        self.record[idx]["attributes"]["convt_time"] = self.convt_time
        self.record[idx]["attributes"]["super_infe"] = self.super_infe