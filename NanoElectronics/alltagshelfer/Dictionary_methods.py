# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 11:19:29 2020

@author: strobe46
"""

a = {'Symm':'yes','cracks':(-100,100,150),'rating':(4,3),'something':'good'}
a.items()
# returns: dict_items([('Symm', 'yes'), ('cracks', (-100, 100, 150)), ('rating', (4, 3)), ('something', 'good')])
list(a.items())
#returns: [('Symm', 'yes'), ('cracks', (-100, 100, 150)), ('rating', (4, 3)), ('something', 'good')]
a.keys()
#returns: dict_keys(['Symm', 'crack', 'rating', 'something'])
list(a.keys())
#returns: ['Symm', 'crack', 'rating', 'something']

'Symm' in a
# retuns: True

'crack' not in a
# returns: False

a.values()
#returns: dict_values(['yes', 'no', (4, 3), 'good'])
list(a.values())
#returns: ['yes', 'no', (4, 3), 'good']

#%% search for items
search = 'good'
for attribute, value in a.items():
    if value == search:
        print(value)
#%%
search = 'cracks'
for attribute, value in a.items():
    if search in attribute:
        print(value)
#%% delete item 
thisdict =	{  "brand": "Ford",  "model": "Mustang",  "year": 1964}
thisdict.pop("model")
print(thisdict) 



