# -*- coding: utf-8 -*-
"""
Created on Fri Apr 06 19:02:13 2018

@author: bmg18
"""




def common_member(a, b):
	a_set = set(a)
	b_set = set(b)
	if (a_set & b_set):
		return list(a_set & b_set)
	else:
	    return []