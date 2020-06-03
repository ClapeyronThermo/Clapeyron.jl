import json

data={"Polymer_length" : {
    "PMMA": 930,
    "PS":961,
    "PB":490,
    "EtOH": 1,
    "Benzene": 1
},

"Segment_mw" : {
    "PMMA": 86,
    "PMA": 100,
    "PS": 104.1,
    "PB": 54.1,
    "EtOH": 46.07,
    "Benzene": 78.11
},

"c_i_coeff" : {
    "PMMA": 0.0092,
    "PS": 0.0066,
},

"R_k" : {
    "CH:CH":1.1167,
    "CH2":0.6744,
    "C":0.2195,
    "CH3COO":1.9031,
    "CH3":0.9011,
    "ACH":0.5313,
    "ACCH":0.8121,
    "OH": 1.0
},

"a_nm" : {
    "CH:CH":{"CH2": 2520,
             "CH:CH":0,
             "C": 2520,
             "CH3": 2520,
             "CH3COO": 71.23,
             "ACH": 340.7,
             "ACCH": 4102,
             "OH": 636.1},
    "CH2":{"CH2": 0,
             "CH:CH":-200,
             "C": 0,
             "CH3": 0,
             "CH3COO": 232.1,
             "ACH": 61.13,
             "ACCH": 76.50,
             "OH":986.5},
    "C":{"CH2": 0,
             "CH:CH":-200,
             "C": 0,
             "CH3": 0,
             "CH3COO": 232.1,
             "ACH": 61.13,
             "ACCH": 76.50,
             "OH": 636.1},
    "CH3COO":{"CH2": 114.8,
             "CH:CH":269.3,
             "C": 114.8,
             "CH3": 114.8,
             "CH3COO": 0,
             "ACH": 85.84,
             "ACCH": -170,
             "OH": 636.1},
    "CH3":{"CH2": 0,
             "CH:CH":-200,
             "C": 0,
             "CH3": 0,
             "CH3COO": 232.1,
             "ACH": 61.13,
             "ACCH": 76.50,
             "OH":986.5},
    "ACH":{"CH2": -11.12,
             "CH:CH":-94.78,
             "C": -11.12,
             "CH3": -11.12,
             "CH3COO": 5.994,
             "ACH": 0,
             "ACCH": 167.0,
             "OH": 636.1},
    "ACCH":{"CH2": -69.7,
             "CH:CH":-269.7,
             "C": -69.7,
             "CH3": -69.7,
             "CH3COO": 5688,
             "ACH": -146.8,
             "ACCH": 0,
             "OH": 636.1},
    "OH":{"CH2": 156.4,
             "CH:CH":0,
             "C": 156.4,
             "CH3": 156.4,
             "CH3COO": 0,
             "ACH": 89.6,
             "ACCH": 0,
             "OH": 0},
},

"QK" : {
    "CH:CH":0.867,
    "CH2":0.540,
    "C":0.0,
    "CH3COO":1.728,
    "CH3":0.848,
    "ACH":0.400,
    "ACCH":0.384,
    "OH": 1.2
},

"Polymer_groups" : {
    "PB":{
        "CH:CH":1,
        "CH2":2
    },
    "PMMA":{
        "CH3":1,
        "CH2":1,
        "C":1,
        "CH3COO":1
    },
    "PS":{
        "ACH": 5,
        "ACCH":1,
        "CH2":1
    },
    "EtOH":{
        "OH": 1,
        "CH3": 1,
        "CH2": 1
    },
    "Benzene":{
        "ACH": 6,
    },
},

"rho_reference" : {
    "PMMA":1.190,
    "PS":1.05,
    "PB":0.915
},

"rho_reference_temp" : {
    "PMMA":20.0,
    "PS":20.0,
    "PB":25.0
},

"T_G" : {
    "PMMA": 105,
    "PS": 80
},

"expansion_coefficients" : {
    "PMMA":{
        "LEQ_TG":2.60e-4,
        "GEQ_TG":5.60e-4
    },
    "PS":{
        "LEQ_TG":1.7e-4,
        "GEQ_TG":5.1e-4
    },
    "PB":6.99e-4
},

"EPSILON" : { 
       "PMMA": 264.6,
       "PS": 348.2,
       "PB": 288.84
},

"SIGMA" : { 
       "PMMA": 3.553e-10,
       "PS": 4.152e-10,
       "PB": 4.097e-10
},

"SEGMENT" : { 
       "PMMA": 0.0270*100,
       "PS": 0.0205*104.1,
       "PB": 0.0245*54.1
},

"Binary_k" : {
"PMMA": {"PS":0.00497},
"PS": {"PMMA":0.00497}

},

"CHI" : {
       "PMMA": {"PS": 0.25},
       "PS": {"PMMA": 0.25}
}}


with open('all_data.json', 'w') as fp:
    json.dump(data, fp)
