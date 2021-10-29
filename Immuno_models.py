import numpy as np
import matplotlib.pyplot as plt
from io import StringIO
from matplotlib.lines import Line2D
from datetime import datetime, timedelta
import scipy.special as sc
import seaborn as sns
import pickle
import json
from scipy.optimize import curve_fit


#----------------- Functions -----------------

def hamming_distance(chaine1, chaine2):

    return sum(c1 != c2 for c1, c2 in zip(chaine1, chaine2))

def find_complementary_seq(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],min(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def find_complementary_seq_2(sequence, Energy_Matrix):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq = list(sequence)
	list_comp_seq = []
	for i in list_seq:
		pos_i = np.where(np.isin(Alphabet,i))[0][0]
		list_comp_seq.append(Alphabet[np.where(np.isin(M[pos_i],max(M[pos_i])))[0][0]])
	comp_seq = "".join(list_comp_seq)
	return comp_seq

def calculate_energy(Energy_Matrix, seq1, seq2):

	M = Energy_Matrix
	Alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't']
	list_seq1 = list(seq1)
	list_seq2 = list(seq2)
	Energy = 0
	for i in range(9):
		pos_i = np.where(np.isin(Alphabet,list_seq1[i]))[0][0]
		pos_j = np.where(np.isin(Alphabet,list_seq2[i]))[0][0]
		Energy += M[pos_i][pos_j]
	return Energy

def my_linear_func(x, a, b):

    return a + b*x

def my_quadratic_func(x, a, b, c):

    return np.log(a)+np.log(np.sqrt(-b)) + b*(x-c)**2

def binding_affinity(x, a, b):
	
    return 1/(1+((np.exp(a+b*x))))

def my_plot_layout(ax, yscale = 'linear', xscale = 'linear', ticks_labelsize = 24, xlabel = '', ylabel = '', title = '', x_fontsize=24, y_fontsize = 24, t_fontsize = 24):
    ax.tick_params(labelsize = ticks_labelsize)
    ax.set_yscale(yscale)
    ax.set_xscale(xscale)
    ax.set_xlabel(xlabel, fontsize = x_fontsize)
    ax.set_ylabel(ylabel, fontsize = y_fontsize)
    ax.set_title(title, fontsize = t_fontsize)

