## converts tables in ipac form to csv
from astropy.table import Table
for i in range(7):
    fname_in = 'radec_extc/' + str(i) + '.tbl'
    fname_out = 'radec_extc/' + str(i) + '.txt'
    t = Table.read(fname_in, format='ipac')
    t.write(fname_out, format='csv')
