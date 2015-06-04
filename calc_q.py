from pylab import real, imag
from matplotlib.pyplot import title, plot, legend, show, savefig
import csv
from os import chdir
import numpy as np
import pandas
import curve

#chdir("D:\Work\Sonnet")

filename = 'A1 gC param'
datatemp = []
f = open(filename + '.csv')
csv_f = csv.reader(f)
print csv_f.line_num

# for row in csv_f:
#     n = 0
#     if len(row) > 1 and n<1:
#         datatemp.append(row)
#     print n
#     n = n+1
#         print datatemp
        #print row

for row in csv_f:
    if len(row) > 1:
        datatemp.append(row)

print 'length of list', len(datatemp)

keys = datatemp.pop(0)
data = np.array(datatemp).astype(np.float)

# print keys
n_freq = keys.index('Frequency (GHz)')
n_reS11 = keys.index('RE[S11]')
n_imS11 = keys.index('IM[S11]')
n_reS21 = keys.index('RE[S21]')
n_imS21 = keys.index('IM[S21]')

freq = data[:,n_freq]
S11 = data[:,n_reS11]+1j*data[:,n_imS11]
S21 = data[:,n_reS21]+1j*data[:,n_imS21]


for label, y, x in (('S11', S11, freq), ('S21',S21, freq)):
    cur = curve.Curve()
    cur.set_data(pandas.Series(y, index=x))
    #res, fit_curve = cur.fit('lorentz_complex_sam')
    res, fit_curve = cur.fit('lorentz_complex_thibault')
    plot(real(fit_curve.data),imag(fit_curve.data), label='fit Q='+str(fit_curve.params['Q_c'])+';f='+str(fit_curve.params['omega_0']))
    #plot(real(fit_curve.data),imag(fit_curve.data), label='fit Q='+str(fit_curve.params['Q'])+';f='+str(fit_curve.params['x0']))
    #plot(real(cur.data),imag(cur.data), 'o', label=label)
legend()
show()
title(filename)

savefig(filename + '.pdf')
savefig(filename + '.png')






