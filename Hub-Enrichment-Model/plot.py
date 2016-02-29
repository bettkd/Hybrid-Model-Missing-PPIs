import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

from hubEnricher import Enricher

data = "Data/Allergy_and_Asthma.txt"
enricher = Enricher(data)
enricher.readPPIsFromFile()
ppis = enricher.ppisdict

# data   
x = np.array(ppis.values())

from sklearn import preprocessing
x = preprocessing.scale(x)

'''
from scipy import stats
x = stats.zscore(x)
'''


mu = np.mean(x) # mean of distribution
sigma = np.std(x) # standard deviation of distribution



num_bins = 285
# the histogram of the data
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='green', alpha=0.75)
# add a 'best fit' line
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.xlabel('Standardized PPI Scores')
plt.ylabel('Frequency of PPI')
plt.title(r'Standardized distribution for PPIs related to Asthma and Allergy')

'''
# Uncomment to plot for distribution
plt.title(r'Distribution of PPI for Asthma and Allergy: $\mu=%.2f$, $\sigma=%.2f$'%(mu, sigma))

'''

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()