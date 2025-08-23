import re

#long_pw = 'PPARA :+: FABP1, FABP1 :+: GPX4, GPX4 :-:  Ferroptosis, FABP1 :-: ACSL4 , ACSL4  :+:  Ferroptosis'
long_pw = 'PPARA :+: FABP1, FABP1 :+: GPX4, GPX4 :-:  Ferroptosis'

interactions = long_pw.split(',')
re.split(r':-:|:\+:', interactions[0])

prev_target = None
all_pws = []
current_pw = []
for reaction in interactions:
    if '+' in reaction:
        source, target = [x.replace(' ', '') for x in re.split(r':\+:', reaction)]
        interaction_type = 1
    if '-' in reaction:
        source, target = [x.replace(' ', '') for x in re.split(r':-:', reaction)]
        interaction_type = -1

    reaction_dict = {
        'source': source,
        'target': target,
        'interaction_type': interaction_type
    }

    if current_pw == []:
        prev_target = target
    if source == prev_target:
        current_pw.append(reaction_dict)
    else:
        if current_pw:
            all_pws.append(current_pw)
        current_pw = [reaction_dict]
    prev_target = target

all_pws.append(current_pw)
