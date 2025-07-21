# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:49:00 2020

@author: devr0035
"""

def supply(self):
    # use bed file
    # zb_file = self.p['bed_file']
    # # get current bathymetry
    zb = self.get_var('zb')
    #now add supply
    zb[:,0:20] += 1
    #zb += 0.0001

    
    # # combine current bathymetry with bed file to compensate waterline interations
    # zb_corr = zb_now
   
    # if zb_file.ndim==1:
    #     zb_corr[0,(zb_now[0,:]-zb_file)<0]=zb_file[(zb_now[0,:]-zb_file)<0]
    # else:
    #     zb_corr[:,(zb_now[0,:]-zb_file)<0]=zb_file[(zb_now[0,:]-zb_file)<0]
   
    # # set net variable
    #print('zb')
    self.set_var('zb',zb)
    
    